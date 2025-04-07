function FParPLScheck2(lambda)

%This code compares different lambda values in arPLS.

%Input
%lambda: the first lambda value used in arPLS. +1, +2, and +3 are also used.

%Output Figure files containing graphs from four different lambda values (a
%value indicated in argument, +1, +2, +3) will be saved.

close all;
functionName = "FParPLScheck2-";

% Parameters

p.frameRate = 30; % in Hz. video frame rate in Doric system.
p.samplingRate = 120; % in Hz. signal sampling rate in Doric system.
p.lambda = lambda; % the first lambda value used in arPLS. +1, +2, and +3 are also used.
p.remove = 15; % in sec. remove first seconds of data to avoid high values at the beginning

% Interactive file selection
disp('Select the input MAT file containing inputInfo:');
[inputFile, inputPath] = uigetfile('*.mat', 'Select Input MAT File');
if isequal(inputFile, 0)
    error('No file selected. Exiting...');
end
inputFile = fullfile(inputPath, inputFile);

% Load information from input file
load(inputFile, 'inputInfo');

% Interactive selection for .doric and .avi files
disp('Select the corresponding .doric files:');
[sigFiles, sigPath] = uigetfile('*.doric', 'Select Signal Files', 'MultiSelect', 'on');
if isequal(sigFiles, 0)
    error('No .doric files selected. Exiting...');
end
if ischar(sigFiles)
    sigFiles = {sigFiles}; % Ensure sigFiles is a cell array
end
sigFiles = fullfile(sigPath, sigFiles);

disp('Select the corresponding .avi files:');
[videoFiles, videoPath] = uigetfile('*.avi', 'Select Video Files', 'MultiSelect', 'on');
if isequal(videoFiles, 0)
    error('No .avi files selected. Exiting...');
end
if ischar(videoFiles)
    videoFiles = {videoFiles}; % Ensure videoFiles is a cell array
end
videoFiles = fullfile(videoPath, videoFiles);

% Check file correspondence
if length(sigFiles) ~= length(videoFiles)
    error('The number of .doric files and .avi files must match.');
end

% Check column "fix" exists and, if so, get information
fixColumn = find(inputInfo.Properties.VariableNames == "fix" | inputInfo.Properties.VariableNames == "Fix");
if length(fixColumn) > 1
    error('There are multiple column "fix" in inputInfo.');
elseif ~isempty(fixColumn)
    fix = inputInfo.fix;
else
    fix = cell(size(inputInfo, 1), 1);
end

% Extract input file name if inputFile contains paths
[~, inputFileName, ~] = fileparts(inputFile);

% Make a folder to save figure but stop if the folder for image has already existed
folderName = strcat(functionName, inputFileName, '-lambda', num2str(p.lambda), '-', string(datetime('now', 'Format', 'yyyy-MM-dd-HH-mm')));
if isfolder(folderName)
    error('Image folder already exists.');
end
mkdir(folderName);

% Run the main function
for ii = 1:length(sigFiles)
    mainfunc(sigFiles{ii}, fix{ii}, p, functionName, folderName);
end
commandwindow;
end

function mainfunc(sigFile, fix, p, functionName, folderName)
    %extract data from .doric file
    Data = ExtractDataAcquisition(sigFile);
    
    %get signal, reference and their time stamps
    [oriSig{1}, oriSig{2}, sigTS] = getFPSignal(Data);
    
    %fix abnormal values indicated by "fix" column in inputInfo by
    %interpolation.
    [oriSig{1}, oriSig{2}, sigTS, nFirstSampleRemoved] = fixAbnormal(oriSig{1}, oriSig{2}, sigTS, fix);
    
    %remove the first p.remove seconds of data at the beggining.
    tmp = p.remove*p.samplingRate - nFirstSampleRemoved;
    oriSig{1}(1:tmp) = [];
    oriSig{2}(1:tmp) = [];
    sigTS(1:tmp) = [];
    
    %get file name to show in figure
    tmp = strfind(sigFile,'\');
    if ~isempty(tmp)
        sigFileName = char(sigFile);
        sigFileName = sigFileName(tmp(end)+1:end);
    else
        sigFileName = sigFile;
    end
    
    %generate figure
    f = figure();
    f.WindowState = 'maximized';
    t = tiledlayout(4,4);
    title(t, sigFileName, 'Interpreter', 'none');
    xlabel(t,'Sec')
    
    sigName = ["Signal vs Baseline","Signal - Baseline","Reference vs Baseline","Reference - Baseline"];
    for ii = 1:4
        tmpL = p.lambda + ii - 1;
        sig{1} = arPLSplus(oriSig{1},10^tmpL,0.001);%baseline of Signal
        sig{2} = oriSig{1} - sig{1};%signal - baseline
        sig{3} = arPLSplus(oriSig{2},10^tmpL,0.001);%baseline of reference
        sig{4} = oriSig{2} - sig{3};%reference - baseline
        for jj = 1:2
            nexttile(4*ii+2*jj-5)
            plot(sigTS,oriSig{jj},"color", [0.5 0.5 0.5]);
            hold on
            plot(sigTS,sig{2*jj-1},"b");
            xlim('tight');
            title(sigName(2*jj-1));
            nexttile(4*ii+2*jj-4)
            plot(sigTS,sig{2*jj},"b");
            xlim('tight');
            title(sigName(2*jj));
        end
    end
    savefig(strcat(folderName,'\',functionName,sigFileName, '_',num2str(p.lambda),'.fig'))
    end
    
    function [outputArg] = ExtractDataAcquisition(filename)
    %filename is the full name of the .doric file where we want to extract all
    %the data from the DataAcquisition
    %
    %OutputArg is a structure with all the Data contained in the .doric file
    
    if ~contains(filename,'.doric')
        filename = [filename '.doric'];
    end
    
    DataAcquisition = h5info(filename,'/DataAcquisition');
    
    
    %Recursive function to go find and extract all the data
        function [Dataset] = getalldata(H5)
            
            if isempty(H5.Groups)
                if ~isempty(H5.Datasets)
                    
                    Dataset_tmp = [];
                    for k = 1:length(H5.Datasets)
                        Name = [H5.Name '/' H5.Datasets(k).Name];
                        Data = h5read(filename,Name);
                        
                        Dataset_tmp(k).Name = H5.Datasets(k).Name;
                        Dataset_tmp(k).Data = Data;
                        Dataset_tmp(k).DataInfo = H5.Datasets(k).Attributes;
                    end
                    
                    Dataset.Name = strrep(H5.Name(2:end),'/','_');
                    Dataset.Data = Dataset_tmp;
                end
            else
                Dataset = [];
                for k=1:length(H5.Groups)
                    Dataset = [Dataset getalldata(H5.Groups(k))];
                end
            end
            
        end
    
    
       outputArg = getalldata(DataAcquisition);
    
    end
    
    function [sig, ref, sigTS] = getFPSignal(Data)
    
    %Check the 1st level of Data, and find 'AIN01xAOUT01-LockIn' and 'AIN01xAOUT02-LockIn'
    dataInd1 = [];
    dataInd2 = [];
    for ii = 1:length(Data)
        dataName = Data(ii).Name;
        if contains(dataName,'AIN01xAOUT01-LockIn')
            dataInd1 = [dataInd1 ii];
        elseif contains(dataName,'AIN01xAOUT02-LockIn')
            dataInd2 = [dataInd2 ii];
        end
    end
    if isempty(dataInd1)
        error('No AIN01xAOUT01-LockIn data in the .doric file.')
    elseif length(dataInd1) > 1
        error('There are multiple AIN01xAOUT01-LockIn data in the .doric file.')
    end
    if isempty(dataInd2)
        error('No AIN01xAOUT02-LockIn data in the .doric file.')
    elseif length(dataInd2) > 1
        error('There are multiple AIN01xAOUT02-LockIn data in the .doric file.')
    end
    
    %Check the 2nd level of Data, and find 'values' and 'Time'.
    dataInd11 = [];
    dataInd12 = [];
    dataInd21 = [];
    dataInd22 = [];
    for ii = 1:length(Data(dataInd1).Data)
        dataName = Data(dataInd1).Data(ii).Name;
        if contains(dataName,'Values') 
            dataInd11 = [dataInd11 ii];
            ref = Data(dataInd1).Data(dataInd11).Data;%reference data
        elseif contains(dataName,'Time') 
            dataInd12 = [dataInd12 ii];
            refTS = Data(dataInd1).Data(dataInd12).Data;%Time stamp of each sample of cameraState
        end
    end
    if isempty(dataInd11)|isempty(dataInd12)
        error('No AIN01xAOUT01-LockIn data in the .doric file.')
    elseif length(dataInd11) > 1|length(dataInd12) > 1
        error('There are multiple AIN01xAOUT01-LockIn data in the .doric file.')
    end
    for ii = 1:length(Data(dataInd2).Data)
        dataName = Data(dataInd2).Data(ii).Name;
        if contains(dataName,'Values') 
            dataInd21 = [dataInd21 ii];
            sig = Data(dataInd2).Data(dataInd21).Data;%main signal data
        elseif contains(dataName,'Time') 
            dataInd22 = [dataInd22 ii];
            sigTS = Data(dataInd2).Data(dataInd22).Data;%Time stamp of each sample of cameraState
        end
    end
    if isempty(dataInd11)|isempty(dataInd12)
        error('No AIN01xAOUT01-LockIn data in the .doric file.')
    elseif length(dataInd11) > 1|length(dataInd12) > 1
        error('There are multiple AIN01xAOUT01-LockIn data in the .doric file.')
    end
    if isempty(dataInd21)|isempty(dataInd22)
        error('No AIN01xAOUT02-LockIn data in the .doric file.')
    elseif length(dataInd11) > 1|length(dataInd12) > 1
        error('There are multiple AIN01xAOUT02-LockIn data in the .doric file.')
    end
    if sigTS ~= refTS
        error('time stamps of signal and reference data do not match.')
    end
    end
    
    %This function check the format and name of data files. It gets animalIDs
    %at the end of file name and check the three files corresponding to each
    %other have the same animalIDs in their file names.
    function animalIDs = checkFileFormat2(sigFile, videoFile)
    animalIDs = nan(size(sigFile,1),2);
    for ii = 1:size(sigFile, 1)
        tmpSigFile = char(sigFile(ii));
        tmpVideoFile = char(videoFile(ii));
    
        if ~isfile(tmpSigFile) || ~isfile(tmpVideoFile)
            error('A specified file does not exist')
        end
    
        if ~strcmp(tmpSigFile(end-5:end),'.doric') || ~strcmp(tmpVideoFile(end-3:end),'.avi')
            error('Check file type in input file. ');
        end
        
        tmpSpaceInd = strfind(tmpSigFile,' ');
        tmpUnderlineInd = strfind(tmpSigFile,'_');
        if isempty(tmpSpaceInd) & ~isempty(tmpUnderlineInd)
            tmpAnimalIDs = tmpSigFile(tmpUnderlineInd(end)+1:end-6);
        elseif ~isempty(tmpSpaceInd) & isempty(tmpUnderlineInd)
            tmpAnimalIDs = tmpSigFile(tmpSpaceInd(end)+1:end-6);
        elseif ~isempty(tmpSpaceInd) & ~isempty(tmpUnderlineInd)
            error('Check animalID part of file names.')
        else
            tmpAnimalIDs = tmpSigFile(max(tmpSpaceInd(end),tmpUnderlineInd(end))+1:end-6);
        end
    
        if ~contains(tmpVideoFile,tmpAnimalIDs)
            error('Check file correspondece between columns in input file.');
        end
    
        tmpHyphenInd = strfind(tmpAnimalIDs, '-');
        if ~isempty(tmpHyphenInd)
            animalIDs(ii,:) = [str2double(tmpAnimalIDs(1:tmpHyphenInd-1)) str2double(tmpAnimalIDs(tmpHyphenInd+1:end))];
        else
            animalIDs(ii,1) = str2double(tmpAnimalIDs);
        end
    end
    end
    
    %fix abnormal values indicated by "fix" column in inputInfo by
    %interpolation.
    function [sig, ref, sigTS, nFirstSampleRemoved] = fixAbnormal(sig, ref, sigTS, fix)
    nFirstSampleRemoved = 0;
    if ~isempty(fix)
        %when signal value is abnormal in the first sample, remove the first
        %consecutive abnormal samples.
        if fix(1) == 1
            tmpInd = find(diff(fix)>1,1,"first");
            if isempty(tmpInd)
                removeInd = fix(end);
            else
                removeInd = fix(tmpInd);
            end
            sig(1:removeInd) = [];
            ref(1:removeInd) = [];        
            sigTS(1:removeInd) = [];
            fix(1:removeInd) = [];
            fix = fix - removeInd;
            nFirstSampleRemoved = tmpInd;
        end
        %when signal value is abnormal in the last sample, remove the last
        %consecutive abnormal samples.
        if ~isempty(fix) & fix(end) == length(sig)
            tmpInd = find(diff(fix)>1,1,"last");
            if isempty(tmpInd)
                removeInd = fix(1);
            else
                removeInd = fix(tmpInd+1);
            end
            sig(removeInd:end) = [];
            ref(removeInd:end) = [];
            sigTS(removeInd:end) = [];
            fix(fix > length(sig)) = [];
        end
        tmpSig1 = sig;
        tmpSig2 =  ref;
        tmpSigTS = sigTS;
        tmpSig1(fix) = [];
        tmpSig2(fix) = [];
        tmpSigTS(fix) = [];
        sig = interp1(tmpSigTS, tmpSig1, sigTS);
        ref = interp1(tmpSigTS, tmpSig2, sigTS);
    end
    end