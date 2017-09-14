function fullAutoFiducializing(inputfilename,inputfiledir)
% this function will replace the 'processAcq' function in the main loop of myProcessing Script

global myScriptData TS XXX


%%%% set up stuff for testing %%%   stuff that will be removed later on
XXX.RunToStartWith = '0137';
XXX.pathToProcessed = '/usr/sci/cibc/Maprodxn/InSitu/17-06-30/Data/Processed'; 
XXX.pathToPreprocessed = '/usr/sci/cibc/Maprodxn/InSitu/17-06-30/Data/Preprocessed';

XXX.FiducialTypes = [2 4 5 7 6];   % the fiducial to process
XXX.nTimeFramesOfOneSearch = 3000;



%%%%%%%%%%%%%%%% here aktuall stuff starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% if XXX not set up, set it up now
if ~isfield(XXX,'firstTime')
    XXX.firstTime = 1;
    setUpStuff
end


%%%% load file into TS, 
curTSidx = loadAndPreprocessFiles(inputfilename,inputfiledir);
myScriptData.curTSidx = curTSidx;

%%%% get RMS and put in TS
RMS=rms(TS{curTSidx}.potvals(XXX.leadsOfAllGroups,:),1);
RMS=RMS-min(RMS);





%%%% get the search areas in current file (the start/end idx in ts to search for beats
nTF = size(TS{curTSidx}.potvals,2);  % number of timeframes in in current TS
searchAreas ={};
startIdx = 1;
endIdx = XXX.nTimeFramesPerSearch;

while endIdx < nTF
    searchAreas{end+1} = [startIdx, endIdx];

    startIdx = startIdx + XXX.nTimeFramesPerSearch;
    endIdx = endIdx + XXX.nTimeFramesPerSearch;
end
endIdx = nTF;
searchAreas{end+1} = [startIdx, endIdx];

% searchAreas ={[1, nTimeFramesPersSeach],[nTimeFramesPersSeach+1, 2*nTimeFramesPersSeach], ... ,[ (n-1)*nTimeFramesPersSeach, nTF] },    the intervals were beats are searched for.
% searchAreas are the indeces in the current potvals that are searched


%%%% main loop: process each searchArea in potvals to get allBeats and allFids
allFids = {};
allBeats = {};

disp([inputfilename '---------------'])
for searchAreaIdx = 1:length(searchAreas)
    
    %%%% get the limmited potvals of that area, 
    searchStartIdx = searchAreas{searchAreaIdx}(1);
    searchEndIdx = searchAreas{searchAreaIdx}(2);
    searchAreaTimeFrames = searchStartIdx:searchEndIdx;
    limPotvalsOfSearchArea = TS{curTSidx}.potvals(XXX.leadsToAutoprocess,searchAreaTimeFrames);
    RMSofSearchArea = RMS(searchAreaTimeFrames);
    
    
    %%%% get templates from userdata, if it is first time 
    if XXX.firstTime % if it is the first file and first area frame to autoprocess
        [XXX.beatKernel,XXX.fidKernels,XXX.locFidValues] = getTemplates(RMS,TS{curTSidx}.potvals); 
        XXX.firstTime = 0;
    end

    
    %%%% now fiducialize searchArea: get fidsOfSearchArea and beatsOfSearchArea based on beatKernel,fidKernels,locFidValues
    [fidsOfSearchArea, beatsOfSearchArea] = fiducializeSearchArea(limPotvalsOfSearchArea,RMSofSearchArea,XXX.beatKernel,XXX.fidKernels,XXX.locFidValues);
 
    
    %%%%% update templates based on lastFoundFids and lastFoundBeats
    if length(beatsOfSearchArea) > XXX.numBeatsToAverageOver  % if enough beats in searchArea (e.g. if not only small searchArea at end of file)
        lastFoundFids = fidsOfSearchArea(end+1-XXX.numBeatsToAverageOver : end);
        lastFoundBeats = beatsOfSearchArea(end+1-XXX.numBeatsToAverageOver : end);
        
        [XXX.beatKernel,XXX.fidKernels,XXX.locFidValues] = updateTemplates(limPotvalsOfSearchArea,RMSofSearchArea, lastFoundFids, lastFoundBeats);
    end
    
    
    
    %%%% put fidsOfSearchArea and beatsOfSearchArea from 'search area frame' in 'complete Run frame' by adding searchStartIdx
    for beatNumber = 1:length(fidsOfSearchArea)
        for fidNumber = 1:length(fidsOfSearchArea{beatNumber})
            fidsOfSearchArea{beatNumber}(fidNumber).value = fidsOfSearchArea{beatNumber}(fidNumber).value + searchStartIdx - 1;
        end
        beatsOfSearchArea{beatNumber} = beatsOfSearchArea{beatNumber} + searchStartIdx - 1;
    end
    
    
    %%%% append the fidsOfSearchArea and beatsOfSearchArea to allFids and allBeats
    allBeats = [allBeats beatsOfSearchArea];
    allFids = [allFids fidsOfSearchArea];

end



%%%%% now the file is fiducialized..   Next step: process each beat and save the files

for beatNumber=1:length(allBeats)
    processBeat(beatNumber,allFids,allBeats)
end









%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [beatKernel,fidKernels,locFidValues] = getTemplates(fullRMS, fullPotvals)
% get the firstTemplates   % this function will later compare many templates to see which one works best

global  myScriptData XXX

%%%% set up stuff
fidsKernelLength = myScriptData.FIDSKERNELLENGTH;
totalKernelLength = 2 * fidsKernelLength +1;   % the length of a kernel
fidTypes=XXX.fidTypes;   % oder here is important: start of a wave must be imediatly followed by end of same wave. otherwise FidsToEvents failes.


%%%% load the tsTempl and the tsPreproc to get the first templates from it
pathToProcessed = XXX.pathToProcessed;
fullPathProcFile = fullfile(pathToProcessed,['Run' XXX.RunToStartWith '-ns.mat']);   % todo: this only works if you start with RunToStartWith!! otherwise wrong template for wrong file


load(fullPathProcFile)
tsTempl = ts;
clear ts






%%%% get beatKernel
bsk = tsTempl.selframes(1);
bek = tsTempl.selframes(2);

beatKernel = fullRMS(bsk:bek);


%%%% get oriFids and remove everything but global fids from it
oriFids=tsTempl.fids;
toBeCleared=[];
for p=1:length(oriFids)
    if length(oriFids(p).value)~=1   % if not single global value
        toBeCleared=[toBeCleared p];
    end
end
oriFids(toBeCleared)=[];


%%%% get locFidValues in the "beat frame" from oriFids
locFidValues = [];
for fidType = fidTypes
    locFidValue = round(oriFids([oriFids.type]==fidType).value);
    locFidValues(end+1) = locFidValue;
end
% get globFidsValues in the "potvals frame" as well
globFidsValues = locFidValues + bsk-1;



%%%% get the first fidKernels based on the user fiducialized beat
fsk=globFidsValues-fidsKernelLength;   % fiducial start kernel,  the index in potvals where the kernel for fiducials starts
fek=globFidsValues+fidsKernelLength;   % analog to fsk, but 'end'

nFids=length(fidTypes);
limPotvals = fullPotvals(XXX.leadsToAutoprocess,:);   % only get kernels for leadsToAutoprocess
nLeads=size(limPotvals,1);
fidKernels = zeros(nLeads,totalKernelLength, nFids);
for fidNumber = 1:nFids
    fidKernels(:,:,fidNumber) = limPotvals(:,fsk(fidNumber):fek(fidNumber));
end
% fidKernels is now nLeads x nTimeFramesOfKernel x nFids array containing all kernels for each lead for each fiducial
% example: kernel(3,:,5)  is kernel for the 3rd lead and the 5th fiducial (in fidTypes)



function [beatKernel,fidKernels,locFidValues] = updateTemplates(limPotvals,RMSofSearchArea, lastFoundFids, lastFoundBeats)
% returns: 
% - fidKernels, a nLeads x nTimeFramesOfKernel x nFids array containing all kernels for each lead for each fiducial
%   example: kernel(3,:,5)  is kernel for the 3rd lead and the 5th fiducial (in fidTypes)
% - beatKernel:  the new kernel to find beat in RMS
% - locFidValues: the values where to search for fids in local 'beat frame' of beatKernel

% inputs:
% - lastFoundFids: the subset of the last nBeats2avrg beats in allFids:    allFids(currentBeat-mBeats2avrg : currentBeat)
% - lastFoundBeats:  the subset of .beats of the last nBeat2avrg:  beats(currentBeat-mBeats2avrg : currentBeat)
% - limPotvals:   limmitted potvals, limPotvals = FullPotvals(framesOfSearchArea, leadsToAutoprocess)


global myScriptData XXX



%%%% set up some stuff
nBeats2avrg = length(lastFoundFids);  % how many beats to average over?
fidsKernelLength = myScriptData.FIDSKERNELLENGTH;
totalKernelLength = 2 * fidsKernelLength +1;   % the length of a kernel
nFids =length(lastFoundFids{1}) / 2;   % how many fids? divide by 2 because there are local and global fids
nLeads =size(limPotvals,1); % how many leads in potvals? 

prevBeatLength = 1 + lastFoundBeats{1}(2) - lastFoundBeats{1}(1);   % the length of previous beats



%%%% first get the new local fids to be found as average of the last nBeats2average beats
allLocalFids = zeros(nBeats2avrg,nFids);   % fids in local "beat frame". 
for beatNum = 1:nBeats2avrg % for each entry in lastFoundFids
    gloFrFidValues = [lastFoundFids{beatNum}(nFids+1:2*nFids).value]; % get fid values of the last already processed beats in global frame
    locFrFidValues = 1 + gloFrFidValues - lastFoundBeats{beatNum}(1);  % put them in local beat frame
    allLocalFids(beatNum,:) = locFrFidValues;   % and store them in allLocFids  , the local frame fid values of the last processed beats
end
% now get the average and round
locFidValues = round(mean(allLocalFids,1));  % the new averaged local fid values in 'old beatKernel frame'


%%%% now average the potential values of the last nBeats2average beats for each lead
allBeatsToAvrg = zeros(nLeads,prevBeatLength,nBeats2avrg);
for beatNum = 1:nBeats2avrg
    timeFramesOfBeat = lastFoundBeats{beatNum}(1):lastFoundBeats{beatNum}(2);
    allBeatsToAvrg(:,:,beatNum) = limPotvals(:,timeFramesOfBeat);
end
% now average over the beats
avrgdPotvalsOfBeat = mean(allBeatsToAvrg, 3);  % averaged individual leads



%%%% now that we have avrgdPtovalsOfBeat and locAvrgdFids, get the new fidKernels
fsk=locFidValues-fidsKernelLength;   % fiducial start kernel,  the index in avrgdPotvalsOfBeat where the kernel for fiducials starts
fek=locFidValues+fidsKernelLength;   % analog to fsk, but 'end'



fidKernels = zeros(nLeads,totalKernelLength, nFids);
for fidNumber = 1:nFids
    fidKernels(:,:,fidNumber) = avrgdPotvalsOfBeat(:,fsk(fidNumber):fek(fidNumber));
end


%%%% get the new beatKernel %%%%%%%%%%

%%%% first get the new beatStart/EndPoints (length of beat changes..)
qrsStartVal = locFidValues(XXX.fidTypes == 2);
tEndVal = locFidValues(XXX.fidTypes == 7);



beatStartShift = qrsStartVal - XXX.distanceFromFid;
beatEndShift = (tEndVal + XXX.distanceFromFid) - prevBeatLength;
newBeatLength = prevBeatLength - beatStartShift + beatEndShift;

%%%%% average RMS over the last nBeats2avrg
allRMSbeats = zeros(nBeats2avrg,newBeatLength);
for beatNum = 1:nBeats2avrg
    newStartIdx = lastFoundBeats{beatNum}(1) + beatStartShift;   % to do: make sure these values dont get out of RMS range during last loop
    newEndIdx = lastFoundBeats{beatNum}(2) + beatEndShift;
    allRMSbeats(beatNum,:) = RMSofSearchArea(newStartIdx:newEndIdx);
end

beatKernel = mean(allRMSbeats,1);

%%%% put the locFidValues from 'old beatKernel frame' in 'new beatKernel frame'
locFidValues = locFidValues - beatStartShift;













function setUpStuff
global XXX myScriptData

XXX.numBeatsToAverageOver = 5;  % how many beats to average over to get new beat template?
XXX.nTimeFramesPerSearch = 10000; % update templates after every nTimeFramesPersSeach time frames

XXX.fidTypes =[2 4 5 7 6]; % the fids to be found. start of wave must follow end of wave

XXX.distanceFromFid = 30;    % newBeatStartIdx = startOfQRS - distanceFromFid,  newBeatEndIdx = endOfTwave + distanceFromFid,  

XXX.updateTemplates = 0;   % dont update at beginning, 


%%%% get the leadsOfAllGroups and filter out badleads
crg=myScriptData.CURRENTRUNGROUP;
badleads=myScriptData.GBADLEADS{crg};     % the global indices of bad leads
leadsOfAllGroups=[myScriptData.GROUPLEADS{crg}{:}];
leadsOfAllGroups=setdiff(leadsOfAllGroups,badleads);  % signal (where the beat is found) will constitute of those.  got rid of badleads
XXX.leadsOfAllGroups = leadsOfAllGroups;  % whenever an RMS is needed, it will consist only of these leadsOfAllgroups


%%%% set leadsToAutoprocess, the leads to find fiducials for and plot.  Only these leads will be used to compute the global fids
nToBeFiducialised=myScriptData.NTOBEFIDUCIALISED;    % nToBeFiducialised  evenly spread leads from leadsOfAllGroups will be chosen for autoprocessing
idxs=round(linspace(1,length(leadsOfAllGroups),nToBeFiducialised));
%idxs=randi([1,length(leadsOfAllGroups)],1,nToBeFiducialised); % this is wrong, since it may create dublicates
XXX.leadsToAutoprocess=leadsOfAllGroups(idxs);





function index = loadAndPreprocessFiles(inputfilename,inputfiledir)
olddir = pwd;
global myScriptData TS myProcessingData;

%%%%% create cellaray files={full acqfilename, mappingfile, calibration file}, if the latter two are needet & exist    
filename = fullfile(inputfiledir,inputfilename);
files{1} = filename;
isMatFile=0;
if contains(inputfilename,'.mat'), isMatFile=1; end


% load & check mappinfile
mappingfile = myScriptData.RUNGROUPMAPPINGFILE{myScriptData.CURRENTRUNGROUP};
if isempty(mappingfile)
    myScriptData.RUNGROUPMAPPINGFILE{myScriptData.CURRENTRUNGROUP} = '';
elseif ~exist(mappingfile,'file')
    msg=sprintf('The provided .mapping file for the Rungroup %s does not exist.',myScriptData.RUNGROUPNAMES{myScriptData.CURRENTRUNGROUP});
    errordlg(msg);
    error('problem with mappinfile.')
else
    files{end+1}=mappingfile;  
end    

if myScriptData.DO_CALIBRATE == 1 && ~isMatFile     % mat.-files are already calibrated
    if ~isempty(myScriptData.CALIBRATIONFILE)
        if exist(myScriptData.CALIBRATIONFILE,'file')
            files{end+1} = myScriptData.CALIBRATIONFILE;
        end
    end
end
    

%%%%%%% read in the files in TS.  index is index with TS{index}=current ts
%%%%%%% structure
    if isMatFile
        index=ioReadMAT(files{:});
    else
        index = ioReadTS(files{:}); % if ac2 file
    end
    
    
%%%%% make ts.filename only the filename without the path

[~,filename,ext]=fileparts(TS{index}.filename);
TS{index}.filename=[filename ext];
    
    
    
    
    
%%%%%% check if dimensions of potvals are correct, issue error msg if not
if size(TS{index}.potvals,1) < myScriptData.MAXLEAD{myScriptData.CURRENTRUNGROUP}
    errordlg('Maximum lead in settings is greater than number of leads in file');
    cd(olddir);
    error('ERROR');
end
cd(olddir)


%%%%  store the GBADLEADS also in the ts structure (in ts.leadinfo)%%%% 
badleads=myScriptData.GBADLEADS{myScriptData.CURRENTRUNGROUP};
TS{index}.leadinfo(badleads) = 1;

%%%%% do the temporal filter of current file %%%%%%%%%%%%%%%%
if myScriptData.DO_FILTER      % if 'apply temporal filter' is selected
    if 0 %isfield(myScriptData,'FILTER')     % this doesnt work atm, cause buttons for Filtersettings etc have been removed
        myScriptData.FILTERSETTINGS = [];
        for p=1:length(myScriptData.FILTER)
            if strcmp(myScriptData.FILTER(p).label,myScriptData.FILTERNAME)
                myScriptData.FILTERSETTINGS = myScriptData.FILTER(p);
            end
        end
    else
        myScriptData.FILTERSETTINGS.B = [0.03266412226059 0.06320942361376 0.09378788647083 0.10617422096837 0.09378788647083 0.06320942361376 0.03266412226059];
        myScriptData.FILTERSETTINGS.A = 1;
    end
    temporalFilter(index);    % no add audit? shouldnt it be recordet somewhere that this was filtered??? TODO
end


function [FidsOfSearchArea, beatsOfSearchArea] = fiducializeSearchArea(potvalsOfSearchArea,RMSofSearchArea,beatKernel,fidKernels,locFidValues)

%%%%% get paramters from myScriptData
global myScriptData XXX


accuracy=myScriptData.ACCURACY;  % abort condition5
fidsKernelLength=myScriptData.FIDSKERNELLENGTH;  % the kernel indices will be from fidsValue-fidsKernelLength  until fidsValue+fidsKernelLength
kernel_shift=0;       % a "kernel shift", to shift the kernel by kernel_shift   % not used, just here as placeholder..
% reminder: it is kernel_idx=fid_start-fidsKernelLength+kernel_shift:fid_start+fidsKernelLength+kernel_shift
window_width=myScriptData.WINDOW_WIDTH;   % dont search complete beat, but only a window with width window_width,
% ws=bs+loc_fidsValues(fidNumber)-window_width;  
% we=bs+loc_fidsValues(fidNumber)+window_width;
fidTypes = XXX.fidTypes;

%%%%% find the beats

beatsOfSearchArea = findMatches(RMSofSearchArea, beatKernel, accuracy);
nBeats=length(beatsOfSearchArea);


%%%% initialice/preallocate allFids
nFids = size(fidKernels,3);
defaultFid(nFids).type=[];
[FidsOfSearchArea{1:nBeats}]=deal(defaultFid);


%%%%%%%%%%%%% fill FidsOfSearchArea with values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=waitbar(0/nBeats,'Autofiducialicing Beats..');
for beatNumber=1:nBeats %for each beat
    bs=beatsOfSearchArea{beatNumber}(1);  % start of beat

    for fidNumber=1:nFids
        %%%% set up windows
        ws=bs+locFidValues(fidNumber)-window_width;  % dont search complete beat, only around fid
        we=bs+locFidValues(fidNumber)+window_width;
        windows=potvalsOfSearchArea(:,ws:we);

            
        
        %%%% find fids
        [globFid, indivFids, variance] = findFid(windows,fidKernels(:,:,fidNumber),'normal');

        %put fids in potval frame
        indivFids=indivFids+fidsKernelLength-kernel_shift+bs-1+locFidValues(fidNumber)-window_width;  % now  newIndivFids is in "complete potvals" frame.
        globFid=globFid+fidsKernelLength-kernel_shift+bs-1+locFidValues(fidNumber)-window_width;      % put it into "complete potvals" frame


        %%%% put the found newIndivFids in FidsOfSearchArea
        FidsOfSearchArea{beatNumber}(fidNumber).type=fidTypes(fidNumber);
        FidsOfSearchArea{beatNumber}(fidNumber).value=indivFids;
        FidsOfSearchArea{beatNumber}(fidNumber).variance=variance;


        %%%% add the global fid to FidsOfSearchArea
        FidsOfSearchArea{beatNumber}(nFids+fidNumber).type=fidTypes(fidNumber);
        FidsOfSearchArea{beatNumber}(nFids+fidNumber).value=globFid; 
    end
    if isgraphics(h), waitbar(beatNumber/nBeats,h), end
end

if isgraphics(h), delete(h), end


function processBeat(beatNumber,allFids,allBeats)
%index: index to orignial ts obtained just before sigSlice in
%myProcessingScript -> mapping, calibration, temporal filter, badleads already done!
%selframes:  frames for slicing  [start:end]
global TS myScriptData XXX


%%%% slice "complete ts" into beat (in TS{newBeatIdx} )
newBeatIdx=tsNew(1);
beatframes=allBeats{beatNumber}(1):allBeats{beatNumber}(2);  % all time frames of the beat

TS{newBeatIdx}=TS{myScriptData.curTSidx};
TS{newBeatIdx}.potvals=TS{newBeatIdx}.potvals(:,beatframes);
TS{newBeatIdx}.numframes=length(beatframes);
TS{newBeatIdx}.selframes=[beatframes(1),beatframes(end)];
    
%%%% put the new fids in the "local beat frame" and save them in newBeatIdx
fids=allFids{beatNumber};
reference=beatframes(1);
for fidNumber=1:length(fids)
    fids(fidNumber).value=fids(fidNumber).value-reference+1;  % fids now in local frame
end
if isfield(fids,'variance'),  fids=rmfield(fids,'variance'); end  %variance not wanted in the output
TS{newBeatIdx}.fids=fids;

%%%%%% if 'blank bad leads' button is selected,   set all values of the bad leads to 0   
if myScriptData.DO_BLANKBADLEADS == 1
    badleads = tsIsBad(newBeatIdx);
    TS{newBeatIdx}.potvals(badleads,:) = 0;
    tsSetBlank(newBeatIdx,badleads);
    tsAddAudit(newBeatIdx,'|Blanked out bad leads');
end


%%%%  baseline correction
if myScriptData.DO_BASELINE
    sigBaseLine(newBeatIdx,[1,length(beatframes)-myScriptData.BASELINEWIDTH],myScriptData.BASELINEWIDTH);
    % also add the baseline fid to ts.fids
    TS{newBeatIdx}.fids(end+1).type=16;
    TS{newBeatIdx}.fids(end).value=1;
    TS{newBeatIdx}.fids(end+1).type=16;
    TS{newBeatIdx}.fids(end).value=length(beatframes)-myScriptData.BASELINEWIDTH;
    
end


%%%%% do activation and deactivation
% if myScriptData.FIDSAUTOACT == 1, DetectActivation(newBeatIdx); end
% if myScriptData.FIDSAUTOREC == 1, DetectRecovery(newBeatIdx); end





%%%% construct the filename  (add eg '-b10' to filename)
[~,filename,~]=fileparts(TS{myScriptData.curTSidx}.filename);
filename=sprintf('%s-b%d',filename,beatNumber); 


%%%% split TS{newIdx} into numGroups smaller ts in grIndices
splitgroup = [];
for p=1:length(myScriptData.GROUPNAME{myScriptData.CURRENTRUNGROUP})
    if myScriptData.GROUPDONOTPROCESS{myScriptData.CURRENTRUNGROUP}{p} == 0, splitgroup = [splitgroup p]; end
end
% splitgroup is now eg [1 3] if there are 3 groups but the 2 should
% not be processed
channels=myScriptData.GROUPLEADS{myScriptData.CURRENTRUNGROUP}(splitgroup);
grIndices = mytsSplitTS(newBeatIdx, channels);    
% update the filenames (add '-groupextension' to filename)
tsDeal(grIndices,'filename',ioUpdateFilename('.mat',filename,myScriptData.GROUPEXTENSION{myScriptData.CURRENTRUNGROUP}(splitgroup))); 
tsClear(newBeatIdx);


%%%% save the new ts structures using ioWriteTS
olddir = cd(myScriptData.MATODIR);
ioWriteTS(grIndices,'noprompt','oworiginal');
cd(olddir);


%%%% do integral maps and save them  
if myScriptData.DO_INTEGRALMAPS == 1
    if myScriptData.DO_DETECT == 0
        msg=sprintf('Need fiducials (at least QRS wave or T wave) to do integral maps for %s.', filename);
        errordlg(msg)
        error('Need fiducials to do integral maps');
    end
    mapindices = fidsIntAll(grIndices);
    if length(splitgroup)~=length(mapindices)
        msg=sprintf('Fiducials (QRS wave or T wave) necessary to do integral maps. However, for %s there are no fiducials for all groups.',filename);
        errordlg(msg)
        error('No fiducials for integralmaps.')
    end

    olddir = cd(myScriptData.MATODIR); 
    fnames=ioUpdateFilename('.mat',filename,myScriptData.GROUPEXTENSION{myScriptData.CURRENTRUNGROUP}(splitgroup),'-itg');

    tsDeal(mapindices,'filename',fnames); 
    tsSet(mapindices,'newfileext','');
    ioWriteTS(mapindices,'noprompt','oworiginal');
    cd(olddir);
    tsClear(mapindices);
end
       
%%%%% Do activation maps   
if myScriptData.DO_ACTIVATIONMAPS == 1
    if myScriptData.DO_DETECT == 0 % 'Detect fiducials must be selected'
        error('Need fiducials to do activation maps');
    end

    %%%% make new ts at TS(mapindices). That new ts is like the old
    %%%% one, but has ts.potvals=[act rec act-rec]
    mapindices = sigActRecMap(grIndices);   


    %%%%  save the 'new act/rec' ts as eg 'Run0009-gr1-ari.mat
    % AND clearTS{mapindex}!
    olddir = cd(myScriptData.MATODIR);
    tsDeal(mapindices,'filename',ioUpdateFilename('.mat',filename,myScriptData.GROUPEXTENSION{myScriptData.CURRENTRUNGROUP}(splitgroup),'-ari')); 
    tsSet(mapindices,'newfileext','');
    ioWriteTS(mapindices,'noprompt','oworiginal');
    cd(olddir);
    tsClear(mapindices);
end

   %%%%% save everything and clear TS
%    saveSettings();          TODO
    tsClear(grIndices);