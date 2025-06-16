%% index back to the whole trace 
ROISegIdx = cell(height(ROI3DWithTraceTable),4);   %create a cell array to store the index of ROI and the corresponding trace
sessionNum = questdlg('Which session of data you want to check:','session check',...
                                   'session1','session2','session3','session1');      %choose the session to check
  
sessionStr = append('registered_trace_',sessionNum);    
wb1 = waitbar(0, 'Indexing the ROI and trace...'); %create a waitbar to show the progress
for i = 1 : height(ROI3DWithTraceTable)
    waitbar(i/height(ROI3DWithTraceTable), wb1); %update the waitbar
    ROISegIdx(i,1) = {i};                       %store the index of ROI
    tempROI = ROI3DWithTraceTable.(sessionStr)(i,1); %get the registered struct of current ROI

    if ~isempty(tempROI.selected_Z)             %if the ROI have registered trace label

        tempROITrace = tempROI.selected_trace;   
        ROISegIdx(i,2) = tempROITrace;         %store the trace of the ROI
    else
        ROISegIdx(i,2) = {0};                  %if the ROI have no registered trace label, store 0
    end
        %assign and calculate the center of each ROI
    tempFP_3D = ROI3DWithTraceTable.FP_3D{i};   %get the footprint of current ROI

    %resize the footprint volume to the same scale as apical dendrite volume (3 times z-step)
    tempFP_3D = imresize3(tempFP_3D,[size(tempFP_3D,1),size(tempFP_3D,2),size(tempFP_3D,3)*3],'nearest');

    [x,y,z] = size(tempFP_3D);                   %get the size of the footprint
    %[roi_x, roi_y,roi_z] = ind2sub([x,y,z],find(tempFP_3D == 1));     %get the index of the footprint, 1 is the positive pixels for this ROI
    %tempCenter = [mean(roi_x), mean(roi_y), mean(roi_z)];              %calculate the center of the footprint
    ROISegIdx(i,3) = table2cell(regionprops3(tempFP_3D,'Centroid'));                                   %store the center of the ROI

end
close(wb1); %close the waitbar

% Remove the ROIs without trace
ROISegIdx(cellfun(@(x) isequal(x, 0), ROISegIdx(:,2)), :) = []; 

%save the ROI
ROISegIdxTable = cell2table(ROISegIdx,'VariableNames',{'ROIIndex','Segmented_trace','Center','TransformedCenter'});