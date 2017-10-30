function file_splitting_unit_test


%%%% get all the folder paths
[pathToUNIT_TESTING_folder,~,~] = fileparts(which('basic_unit_test.m'));   % find path to PFEIFER.m

pathToScriptDataFile = fullfile(pathToUNIT_TESTING_folder,'FILE_SPLITTING_UNIT_TEST','fileSplittingHelperFiles','fileSplittingUnitTestScriptData.mat');
pathToTemplateOutput = fullfile(pathToUNIT_TESTING_folder,'FILE_SPLITTING_UNIT_TEST','templateFolderFileSplitUnitTest');
testOutputDir =  fullfile(pathToUNIT_TESTING_folder,'FILE_SPLITTING_UNIT_TEST','outputFolderFileSplittingUnitTest');
testInputDir = fullfile(pathToUNIT_TESTING_folder,'FILE_SPLITTING_UNIT_TEST','inputFolderFileSplittingUnitTest');


%%%%% start PFEIFER and get figure
PFEIFER
settingsFigure = findobj(allchild(0),'tag','PROCESSINGSCRIPTSETTINGS');

%%%% load the helper files into PFEIFER  (simulate user loading the helper files
try
    executeCallback(settingsFigure,'SCRIPTFILE',pathToScriptDataFile)
catch
    disp('-----UNIT TEST: file_splitting_unit_test --------')
    disp('ERROR: Could not load helper files')
end


%%%%% set input path 
executeCallback(settingsFigure,'ACQDIR',testInputDir)

%%%% open the File Splitting window and get the figure obj
selectSplitfiles % aequivalent to mouse click on FILE SPLITTING
splitFileFigure = findobj(allchild(0),'Name','FILE SPLITTING');


%%%% set some parameters in split file window
executeCallback(splitFileFigure,'SPLITDIR',testOutputDir)
executeCallback(splitFileFigure,'SPLITINTERVAL','3')
executeCallback(splitFileFigure,'SELECTALL','dummy')
executeCallback(splitFileFigure,'CALIBRATE_SPLIT',0)


%%%% execute script -  simulate mouseclick on SPLIT FILES
selectSplitfiles('splitFiles',splitFileFigure)



%%%% compare output with template output 
feedbackMsg = compareFolders(testOutputDir,pathToTemplateOutput);
if ~strcmp(feedbackMsg,'foldersAreEqual')
    disp('-----UNIT TEST: file_splitting_unit_test --------')
    disp('ERROR: The PFEIFER output folder looked different than the template folder. This is the feedbackMsg msg:')
    fprintf(feedbackMsg);
end

%%%% clear outputDir
clearFolder(testOutputDir)

%%%% close matmap
closePFEIFER


disp('SUCCESSFULLY RAN UNIT TEST: file_splitting_unit_test')



function executeCallback(fig,tag,input)
hObject =  findobj(allchild(fig),'Tag',tag);
callbackFkt = hObject.Callback;

switch tag
    case {'SAVESCRIPTDATA','SAVEPROCESSINGDATA'}
        callbackFkt(hObject,input)
    otherwise
        switch hObject.Style
            case {'slider','radiobutton','togglebutton', 'checkbox','listbox', 'popupmenu'}
                try
                    hObject.Value = input;
                catch
                    error('UNITTEST Could not execute callback with tag %s in figure %s.', tag,fig.Name)
                end
            case {'edit'}
                hObject.String = input;
        end
        callbackFkt(hObject,'dummyEventData');   
end



function closePFEIFER
settingsFigure = findobj(allchild(0),'tag','PROCESSINGSCRIPTSETTINGS');
if isempty(settingsFigure), return, end
PFEIFER('CloseFcn',settingsFigure)


function feedbackMsg = compareFolders(folder1, folder2)
feedbackMsg = 'foldersAreEqual';
%%%%% first test if the folders have the same files
metaFolder1 = dir(fullfile(folder1, '*.mat'));
metaFolder2 = dir(fullfile(folder2, '*.mat'));
metaFolder1( [metaFolder1.isdir] ) = [];
metaFolder2( [metaFolder2.isdir] ) = [];
files1={metaFolder1.name};
files2 = {metaFolder2.name};
files_not_in_folder2  = setdiff( files1, files2 );
files_not_in_folder1       = setdiff( files2, files1);

if ~isempty(files_not_in_folder2)
    feedbackMsg = sprintf('There is this files exclusively in folder %s:\n %s\n',folder1,files_not_in_folder2{1});
    return
elseif ~isempty(files_not_in_folder1)
    feedbackMsg = sprintf('There are these files exclusively in folder %s:\n %s\n',folder2,files_not_in_folder1{1});
    return
end

%%%%% if files are the same, load them and compare contents:
for fileIdx = 1:length(files1)
    path1 = fullfile(folder1,files1{fileIdx});
    path2 = fullfile(folder2,files1{fileIdx});
    
    
    fileMeta1 = load(path1);
    fileMeta2 = load(path2);
    
    
    %%%% compare if they have same variables
    
    varNames1=fieldnames(fileMeta1);
    varNames2=fieldnames(fileMeta2);
    if ~isequal(varNames1, varNames2)
        feedbackMsg = sprintf('the files %s don''t have same variables in both folders. \n',files1{fileIdx});
        return
    end
    
    
    
    %%%% compare the ts
    fieldsToIgnore = {'audit','potvals'};
    pv1=fileMeta1.ts.potvals;
    pv2=fileMeta2.ts.potvals;
    
    if ~isempty(find( (pv1-pv2) > 0.001, 1))
        feedbackMsg = sprintf('the files %s don''t have the same potvals ... \n',files1{fileIdx});
    end
    
    
    ts1=rmfield(fileMeta1.ts,fieldsToIgnore);
    ts2=rmfield(fileMeta2.ts,fieldsToIgnore);
    if ~isequal(ts1, ts2)
        feedbackMsg = sprintf('the files %s don''t have the same ts ... \n',files1{fileIdx});
        return
    end
end


function clearFolder(pathToFolder)
curDir = cd(pathToFolder);
delete *.mat
cd(curDir)
