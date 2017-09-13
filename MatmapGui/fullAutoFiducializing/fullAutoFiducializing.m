function fullAutoFiducializingOfFile(inputfilename,inputfiledir)
% this function will replace the 'processAcq' function in the main loop of myProcessing Script

global myScriptData TS XXX


%%%% set up stuff for testing %%%   stuff that will be removed later on
XXX.RunToStartWith = '0137';
XXX.pathToProcessed = '/usr/sci/cibc/Maprodxn/InSitu/17-06-30/Data/Processed'; 

XXX.FiducialTypes = [2 4 5 7 6];   % the fiducial to process
XXX.nTimeFramesOfOneSearch = 10000;


%%%% if XXX not set up, set it up now
if ~isfield(XXX,'settedUp')
    XXX.settedUp = 1;
    setUpStuff
end



%%%% get the template beat for autofiducializing
if XXX.FirstFile % if first time, no templates set yet
    getTemplates    
    XXX.FirstFile = 0;
end



%%%% load file into TS, 
TSidx = loadAndPreprocessFiles(inputfilename,inputfiledir);

%%%% get RMS and put in TS
RMS=rms(TS{TSidx}.potvals(XXX.leadsOfAllGroups,:),1);
RMS=RMS-min(RMS);

RMSidx = tsNew;
TS{RMSidx} = RMS;




%%%% get the search areas in current file (the start/end idx in ts to search for beats
nTFprev = nTFcur;   %number of Timeframes in prev and cur ts
nTFcur = size(TS{curTSindex}.potvals,2);

startIdx = 1;
endIdx = XXX.nTimeFramesPerSearch;
if XXX.nTimeFramesPerSearch > nTFcur
    while endIdx < nTFcur
        searchAreas{end+1} = [startIdx, endIdx];

        startIdx = startIdx + XXX.nTimeFramesPerSearch;
        endIdx = endIdx + XXX.nTimeFramesPerSearch;
    end
    endIdx = XXX.nTimeFramesPerSearch;
    searchAreas{end+1} = [startIdx, endIdx];
else  % if nFrames is less then nTimeFramesPerSearch => if file is to short to be splitted in intervals
    searchAreas{end+1} = [1, nTFcur];
end
% searchAreas ={[1, nTimeFramesPersSeach],[nTimeFramesPersSeach+1, 2*nTimeFramesPersSeach], ... ,[ (n-1)*nTimeFramesPersSeach, nTF] },    the intervals were beats are searched for.
% searchAreas are the indeces in the current potvals that are searched


%%%% main loop: process each searchArea in potvals to get allBeats and allFids
allFids = {};
allBeats = {};
for searchAreaIdx = 1:length(searchAreas)
    
    %%%% update beatKernel,fidKernels,locFidValues  based on the last nBeatsToAverage fiducialized beats of last searchArea
    if XXX.updateTemplates    % if the last interval (possibly from previous file) had less then nBeatsToAverage beats, dont upgrade template
        updateTemplates    % update beatTemplate and fidTemplate
    end


    
    %%%% now fiducialize searchArea: get fidsOfSearchArea and beatsOfSearchArea based on beatKernel,fidKernels,locFidValues
    searchStartIdx = searchAreas{searchAreaIdx}(1);
    searchEndIdx = searchAreas{searchAreaIdx}(2);
    searchAreaTimeFrames = searchStartIdx:searchEndIdx;
    potvalsToInvestigate = potvals(XXX.leadsToAutoprocess,searchAreaTimeFrames);
    [fidsOfSearchArea, beatsOfSearchArea] = fiducializeSearchArea(potvalsToInvestigate,RMS,beatKernel,fidKernels,locFidValues);
    
    
    
    %%%% put fidsOfSearchArea and beatsOfSearchArea from 'search area frame' in 'complete Run frame' by adding searchStartIdx
    for beatNumber = 1:lenght(fidsOfSearchArea)
        for fidNumber = 1:length(fidsOfSearchArea{beatNumber}.value)
            fidsOfSearchArea{beatNumber}(fidNumber).value = fidsOfSearchArea{beatNumber}(fidNumber).value + searchStartIdx - 1;
        end
        beatsOfSearchArea{beatNumber} = beatsOfSearchArea{beatNumber} + searchStartIdx - 1;
    end
    
    
    %%%% append the fidsOfSearchArea and beatsOfSearchArea to allFids and allBeats
    allBeats = [allBeats beatsOfSearchArea];
    allFids = [allFids fidsOfSearchArea];

end



%%%%% now the file is fiducialized..   Next step: process each beat and save the files

for beatNumber=1:length(AUTOPROCESSING.beats)
    processBeat(beatNumber)
end









%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [beatKernel,fidKernels,locFidValues] = getTemplates
% get the firstTemplates   % this function will later compare many templates to see which one works best

global  myScriptData XXX

%%%% set up stuff

totalKernelLength = 2 * myScriptData.fidsKernelLength +1;   % the length of a kernel
fidsTypes=[2 4 5 7 6];   % oder here is important: start of a wave must be imediatly followed by end of same wave. otherwise FidsToEvents failes.


%%%% load the ts to get the first templates from it
pathToProcessed = XXX.pathToProcessed;
fullPathToTemplateFile = [pathToProcessed 'Run' XXX.RunToStartWith '-ns.mat'];
load(fullPathToTemplateFile)


%%%% get beatKernel
bsk = ts.selframes(1);
bek = ts.selframes(2);

RMS=rms(ts.potvals(XXX.leadsOfAllGroups,:),1);
RMS=RMS-min(RMS);

beatKernel = RMS(bsk:bek);


%%%% get oriFids and remove everything but global fids from it
oriFids=ts.fids;
toBeCleared=[];
for p=1:length(oriFids)
    if length(oriFids(p).value)~=1   % if not single global value
        toBeCleared=[toBeCleared p];
    end
end
oriFids(toBeCleared)=[];


%%%% get locFidValues in the "beat frame" from oriFids
for fidType = fidsTypes
    locFidValue = round(oriFids([oriFids.type]==fidType).value);
    locFidValues(end+1) = locFidValue;
end
% get globFidsValues in the "potvals frame" as well
globFidsValues = locFidValues + bsk-1;



%%%% get the first fidKernels based on the user fiducialized beat
fsk=globFidsValues-fidsKernelLength+kernel_shift;   % fiducial start kernel,  the index in potvals where the kernel for fiducials starts
fek=globFidsValues+fidsKernelLength+kernel_shift;   % analog to fsk, but 'end'

nFids=length(fidsTypes);
potvals = ts.potvals(XXX.leadsToAutoprocess,:);   % only get kernels for leadsToAutoprocess
nLeads=size(potvals,1);
fidKernels = zeros(nLeads,totalKernelLength, nFids);
for fidNumber = 1:nFids
    fidKernels(:,:,fidNumber) = potvals(:,fsk(fidNumber):fek(fidNumber));
end
% fidKernels is now nLeads x nTimeFramesOfKernel x nFids array containing all kernels for each lead for each fiducial
% example: kernel(3,:,5)  is kernel for the 3rd lead and the 5th fiducial (in fidsTypes)









function updateTemplates
% get a new beatTemplate (fids and beatkernel) based on the last XXX.NumToAverageOver beats



function setUpStuff
global XXX myScriptData

XXX.beatInterval = 30;    % after how many beats should template be updated? 
XXX.numToAverageOver = 5;  % how many beats to average over to get new beat template?
XXX.nTimeFramesPerSearch = 10000; % update templates after every nTimeFramesPersSeach time frames

XXX.lastProcTimeFrameLastFile =[]; 
XXX.startTimeFrameCurFile = 1;
XXX.endTimeFrameCurFile = 1;

XXX.overlapWithPrev = 0;
XXX.overlapWithNext = 0;

XXX.beatsPrev = {};
XXX.beatsCurr = {};

XXX.beatsPrevStartIdx = 0;
XXX.beatsCurStartIdx = 0;
XXX.beatsCurEndIdx = 0;






%%%% get the leadsOfAllGroups and filter out badleads
crg=myScriptData.CURRENTRUNGROUP;
badleads=myScriptData.GBADLEADS{crg};     % the global indices of bad leads
leadsOfAllGroups=[myScriptData.GROUPLEADS{crg}{:}];
leadsOfAllGroups=setdiff(leadsOfAllGroups,badleads);  % signal (where the beat is found) will constitute of those.  got rid of badleads
XXX.leadsOfAllGroups = leadsOfAllGroupsOfAllGroups;  % whenever an RMS is needed, it will consist only of these leadsOfAllgroups


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


function [FidsOfSearchArea, beatsOfSearchArea] = fiducializeSearchArea(potvals,RMS,beatKernel,fidKernels,locFidValues)

%%%%% get paramters from myScriptData
global myScriptData

accuracy=myScriptData.ACCURACY;  % abort condition5
fidsKernelLength=myScriptData.FIDSKERNELLENGTH;  % the kernel indices will be from fidsValue-fidsKernelLength  until fidsValue+fidsKernelLength
kernel_shift=0;       % a "kernel shift", to shift the kernel by kernel_shift   % not used, just here as placeholder..
% reminder: it is kernel_idx=fid_start-fidsKernelLength+kernel_shift:fid_start+fidsKernelLength+kernel_shift
window_width=myScriptData.WINDOW_WIDTH;   % dont search complete beat, but only a window with width window_width,
% ws=bs+loc_fidsValues(fidNumber)-window_width;  
% we=bs+loc_fidsValues(fidNumber)+window_width;


%%%%% find the beats
beatsOfSearchArea=findMatches(RMS, beatKernel, accuracy);
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
        windows=potvals(:,ws:we);        
        
        %%%% find fids
        [globFid, indivFids, variance] = findFid(windows,fidKernels(:,:,fidNumber),'normal');

        %put fids in potval frame
        indivFids=indivFids+fidsKernelLength-kernel_shift+bs-1+locFidValues(fidNumber)-window_width;  % now  newIndivFids is in "complete potvals" frame.
        globFid=globFid+fidsKernelLength-kernel_shift+bs-1+locFidValues(fidNumber)-window_width;      % put it into "complete potvals" frame


        %%%% put the found newIndivFids in FidsOfSearchArea
        FidsOfSearchArea{beatNumber}(fidNumber).type=fidsTypes(fidNumber);
        FidsOfSearchArea{beatNumber}(fidNumber).value=indivFids;
        FidsOfSearchArea{beatNumber}(fidNumber).variance=variance;


        %%%% add the global fid to FidsOfSearchArea
        FidsOfSearchArea{beatNumber}(nFids+fidNumber).type=fidsTypes(fidNumber);
        FidsOfSearchArea{beatNumber}(nFids+fidNumber).value=globFid; 
    end
    if isgraphics(h), waitbar(beatNumber/nBeats,h), end
end

if isgraphics(h), delete(h), end


function processBeat(beatNumber)
%index: index to orignial ts obtained just before sigSlice in
%myProcessingScript -> mapping, calibration, temporal filter, badleads already done!
%selframes:  frames for slicing  [start:end]
global TS myScriptData XXX


%%%% slice "complete ts" into beat (in TS{newBeatIdx} )
newBeatIdx=tsNew(1);
beatframes=XXX.beats{beatNumber}(1):XXX.beats{beatNumber}(2);  % all time frames of the beat

TS{newBeatIdx}=TS{myScriptData.unslicedDataIndex};
TS{newBeatIdx}.potvals=TS{newBeatIdx}.potvals(:,beatframes);
TS{newBeatIdx}.numframes=length(beatframes);
TS{newBeatIdx}.selframes=[beatframes(1),beatframes(end)];
    
%%%% put the new fids in the "local beat frame" and save them in newBeatIdx
fids=XXX.allFids{beatNumber};
reference=beatframes(1);
for fidNumber=1:length(fids)
    fids(fidNumber).value=fids(fidNumber).value-reference+1;  % fids now in local frame
end
if isfield(fids,'variance'),  fids=rmfield(fids,'variance'); end  %variance not wanted in the output
TS{newBeatIdx}.fids=fids;


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
if myScriptData.FIDSAUTOACT == 1, DetectActivation(newBeatIdx); end
if myScriptData.FIDSAUTOREC == 1, DetectRecovery(newBeatIdx); end





%%%% construct the filename  (add eg '-b10' to filename)
[~,filename,~]=fileparts(TS{myScriptData.unslicedDataIndex}.filename);
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