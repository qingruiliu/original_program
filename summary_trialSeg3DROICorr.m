%% Get the segmented calcium trace from the ROI3DWithTraceTable of each ROI

%select and load the ROICorrIdxTable from the first animal(M3)
msgbox('select the ROICorrIdxTable.mat file for the first animal(M3)');
pause(2)
[name,path] = uigetfile('ROICenterAnd*.mat'); %load the center and transformed center table of M3
cd(path);
load(name);

%select the segmented calcium trace table
disp('---select the ROI3DWithTraceTable.mat file for the first animal(M3)---');
[name,path] = uigetfile('ROI3DWithTraceTable*.mat'); %load the segmented calcium trace table of M3
cd(path);
load(name);

%replace the trace column in ROICorrIdxTable with the segmented trace
isTraceIdx = ROICorrIdxTable.ROIIndex;

for i = 1: length(isTraceIdx)
    tempIdx = isTraceIdx(i);
    ROICorrIdxTable.Trace{i} = ROI3DWithTraceTable.S1_registered(tempIdx,1).selected_trace;
end

%rename and save the variable column in the new table
ROISegTraceTable = renamevars(ROICorrIdxTable,'Trace','Segmented_trace');
save('ROISegTraceTable.mat','ROISegTraceTable','-v7.3');
clearvars -except ROISeg*

%% select and load the ROICorrIdxTable from the second animal(M9)
msgbox('select the ROICorrIdxTable.mat file for the first animal(M9)');
pause(2)
[name,path] = uigetfile('ROICenterAnd*.mat'); %load the center and transformed center table of M3
cd(path);
load(name);

%select the segmented calcium trace table
disp('---select the ROI3DWithTraceTable.mat file for the first animal(M9)---');
[name,path] = uigetfile('ROI3DWithTraceTable*.mat'); %load the segmented calcium trace table of M3
cd(path);
load(name);

%replace the trace column in ROICorrIdxTable with the segmented trace
isTraceIdx = ROICorrIdxTable.ROIIndex;

for i = 1: length(isTraceIdx)
    tempIdx = isTraceIdx(i);
    ROICorrIdxTable.Trace{i} = ROI3DWithTraceTable.S1_registered(tempIdx,1).selected_trace;
end

%rename and save the variable column in the new table
ROISegTraceTable = renamevars(ROICorrIdxTable,'Trace','Segmented_trace');
save('ROISegTraceTable.mat','ROISegTraceTable','-v7.3');
clearvars -except ROISeg*

%% select the tables from 3rd animal
msgbox('select the ROICorrIdxTable.mat file for the first animal(M17)');
pause(2)
[name,path] = uigetfile('ROICenterAnd*.mat'); %load the center and transformed center table of M3
cd(path);
load(name);

%select the segmented calcium trace table
disp('---select the ROI3DWithTraceTable.mat file for the first animal(M17)---');
[name,path] = uigetfile('ROI3DWithTraceTable*.mat'); %load the segmented calcium trace table of M3
cd(path);
load(name);

%replace the trace column in ROICorrIdxTable with the segmented trace
isTraceIdx = ROICorrIdxTable.ROIIndex;

for i = 1: length(isTraceIdx)
    tempIdx = isTraceIdx(i);
    ROICorrIdxTable.Trace{i} = ROI3DWithTraceTable.S1_registered(tempIdx,1).selected_trace;
end

%rename and save the variable column in the new table
ROISegTraceTable = renamevars(ROICorrIdxTable,'Trace','Segmented_trace');
save('ROISegTraceTable.mat','ROISegTraceTable','-v7.3');
clear all

%% start analysis
%ask for the important parameter settings
dlgTitle = 'Set the threshold parameters';
promptPara = {'Set the CC threshold','Set the spatial distance threshold'};
fieldSize = [1 100;1 100];
defInput = {'0.5','20'};
thresholds = inputdlg(promptPara,dlgTitle,fieldSize,defInput);

ccThreshold = str2double(thresholds{1});
spatialThreshold = str2double(thresholds{2});

%% load the ROISegTraceTable from 3 animals
%M3
disp('----------choose the ROISegTraceTable.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M3 = ROISegTraceTable;

%M9
disp('----------choose the ROISegTraceTable.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M9 = ROISegTraceTable;

%M17
disp('----------choose the ROISegTraceTable.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M17 = ROISegTraceTable;


%% next step: contaneate the traces in 100% hit trials 
% (and limit the trial condition in the trial before the hit trial)