%% Based on the cellReg registrated result,assign the 2P mask to the according index

%choose the session currently working with
sessionNum = questdlg('Which session of data are you working with:','session check',...
                                   'session1','session2','session3','session1');      %choose the session to check
  
sessionStr = append('registered_trace_',sessionNum);        %11.5 update, choose the session working with
% load the cellReg result and only keep the information with 2P info
disp('----------choose the cellReg output .mat----------')
[CROut,CRpath] = uigetfile('cellRegistered*.mat');
cd(CRpath)
load(CROut)
cellRegRegisterMap = cell_registered_struct.cell_to_index_map;
nullIdx = cellRegRegisterMap(:,1) == 0;
cellRegRegisterMap(nullIdx,:) = [];
 %important update, cellReg output may neglect the order in the first
 %column
cellRegRegisterMap = sortrows(cellRegRegisterMap,1,'ascend'); 

%% 1st column: 2P mask; later columns: cellpose mask
disp('--------choose the 2P mat file--------')
[refNum,path] = uigetfile('cellRegIn*.mat');
cd(path)
twoPName = extractBefore(refNum,'.');

disp('--------choose the 1st reference plane--------')
[refNum,path] = uigetfile('Z*.mat');
cd(path)
refName1 = extractBefore(refNum,'.');

disp('--------choose the 2nd reference plane--------')
[refNum,path] = uigetfile('Z*.mat');
cd(path)
refName2 = extractBefore(refNum,'.');

disp('--------choose the 3rd reference plane--------')
[refNum,path] = uigetfile('Z*.mat');
cd(path)
refName3 = extractBefore(refNum,'.');

disp('--------choose the 4th reference plane--------')
[refNum,path] = uigetfile('Z*.mat');
cd(path)
refName4 = extractBefore(refNum,'.');

disp('--------choose the 5th reference plane--------')
[refNum,path] = uigetfile('Z*.mat');
cd(path)
refName5 = extractBefore(refNum,'.'); 

%% give the variable of each column's name

registerMapTitle = {twoPName, refName1, refName2, refName3, refName4, refName5};

registerMapTable = array2table(cellRegRegisterMap,'VariableNames',registerMapTitle);

%create a new index map, convert the cellReg index to the 3D-ROI index
disp('----------choose the suite2pTable.mat of current plane----------')
[S2POut,S2Ppath] = uigetfile('suite2p*.mat');
cd(S2Ppath)
load(S2POut)
register3DROITable = registerMapTable;
register3DROITable.(twoPName)(1:size(register3DROITable,1),1) = suite2pTable.ROIindexS2P((1:size(register3DROITable,1)));

%% use the registered map to refer the corresponding 3D ROI
%get the index number (reference plane 1)
refPlane1 = extractAfter(refName1,'Z');                             %get the number of reference plane #1
refIdx1 = find(coordinateTable.planeN == str2double(refPlane1));    %get the number of the reference plane in the coordinateTable

refCoorIdxList = coordinateTable.ROIIdxList{refIdx1};               %3D ROI value
refCoorIdx = registerMapTable.(refName1);                           %index of 3D ROI value
tempOutIdx = zeros(size(refCoorIdx));
nZIdx = refCoorIdx ~= 0;                                            
tempOutIdx(nZIdx) = refCoorIdxList(refCoorIdx(nZIdx));              %give value to the non-zero elment using the non-zero index

register3DROITable.(refName1) = tempOutIdx;

% get the index number (reference plane 2)
refPlane2 = extractAfter(refName2,'Z');                             %get the number of reference plane #1
refIdx2 = find(coordinateTable.planeN == str2double(refPlane2));    %get the number of the reference plane in the coordinateTable

refCoorIdxList = coordinateTable.ROIIdxList{refIdx2};               %3D ROI value
refCoorIdx = registerMapTable.(refName2);                           %index of 3D ROI value
tempOutIdx = zeros(size(refCoorIdx));
nZIdx = refCoorIdx ~= 0;                                            
tempOutIdx(nZIdx) = refCoorIdxList(refCoorIdx(nZIdx));              %give value to the non-zero elment using the non-zero index

register3DROITable.(refName2) = tempOutIdx;

% get the index number (reference plane 3)
refPlane3 = extractAfter(refName3,'Z');                             %get the number of reference plane #1
refIdx3 = find(coordinateTable.planeN == str2double(refPlane3));    %get the number of the reference plane in the coordinateTable

refCoorIdxList = coordinateTable.ROIIdxList{refIdx3};               %3D ROI value
refCoorIdx = registerMapTable.(refName3);                           %index of 3D ROI value
tempOutIdx = zeros(size(refCoorIdx));
nZIdx = refCoorIdx ~= 0;                                            
tempOutIdx(nZIdx) = refCoorIdxList(refCoorIdx(nZIdx));              %give value to the non-zero elment using the non-zero index

register3DROITable.(refName3) = tempOutIdx;

% get the index number (reference plane 4)
refPlane4 = extractAfter(refName4,'Z');                             %get the number of reference plane #1
refIdx4 = find(coordinateTable.planeN == str2double(refPlane4));    %get the number of the reference plane in the coordinateTable

refCoorIdxList = coordinateTable.ROIIdxList{refIdx4};               %3D ROI value
refCoorIdx = registerMapTable.(refName4);                           %index of 3D ROI value
tempOutIdx = zeros(size(refCoorIdx));
nZIdx = refCoorIdx ~= 0;                                            
tempOutIdx(nZIdx) = refCoorIdxList(refCoorIdx(nZIdx));              %give value to the non-zero elment using the non-zero index

register3DROITable.(refName4) = tempOutIdx;

% get the index number (reference plane 5)
refPlane5 = extractAfter(refName5,'Z');                             %get the number of reference plane #1
refIdx5 = find(coordinateTable.planeN == str2double(refPlane5));    %get the number of the reference plane in the coordinateTable

refCoorIdxList = coordinateTable.ROIIdxList{refIdx5};               %3D ROI value
refCoorIdx = registerMapTable.(refName5);                           %index of 3D ROI value
tempOutIdx = zeros(size(refCoorIdx));
nZIdx = refCoorIdx ~= 0;                                            
tempOutIdx(nZIdx) = refCoorIdxList(refCoorIdx(nZIdx));              %give value to the non-zero elment using the non-zero index

register3DROITable.(refName5) = tempOutIdx;

%% assign the DeconvF value to the 3D ROI table(suite2pTable --> registered3DROITable --> ROI3DWithTraceTable)

% get the most frequent non-zero element in register3DROITable

ROIIdxMat = table2array(register3DROITable(:,2:end));   %exclude the first column(2P ROI index)

mostNonZero = zeros(size(ROIIdxMat,1),1);               %create the output array

for i = 1 : size(ROIIdxMat)

    tempRow = ROIIdxMat(i,:);                           %extract the values of the current row
    nonZeroValues = tempRow(tempRow ~= 0);              %get all of the non-zero elements

    if ~isempty(nonZeroValues)
        [uniqueVals, ~, idx] = unique(nonZeroValues);  %get unique value and idx mapping
        freqCount = accumarray(idx,1);                 %count occurrences
        [~,maxIdx] = max(freqCount);                   %get the index of max frequent
        mostNonZero(i) = uniqueVals(maxIdx);           %get the most frequent non-zero value
        
        %assign the 2P information to 3D ROI table
        tempIdx = mostNonZero(i);                      %most frequent non-zero element is the 3D ROI index
        tempTrace = suite2pTable.DeconvF(i);           %get the deconvoluted trace of 2P ROI
            
        depthLabel = extractAfter(twoPName,'t');       %current 2P plane label
        ROI3DWithTraceTable.(sessionStr)(tempIdx).(depthLabel) = tempTrace;   %10.31 session2
    else
        mostNonZero(i) = 0;
    end

end
%% 
registerTableName = append('register3DROITable',depthLabel,'Binary.mat');
save(registerTableName,'register3DROITable','-v7.3')
clearvars -except *Table