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
clearvars -except ROISeg*

%% load the ROISegTraceTable from each animal and prepare for the analysis