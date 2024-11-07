%% index back to the whole trace 
ROICorrIdx = cell(height(ROI3DWithTraceTable),2);   %create a cell array to store the index of ROI and the corresponding trace
for i = 1 : height(ROI3DWithTraceTable)
    ROICorrIdx(i,1) = {i};                          %store the index of ROI
    tempROI = ROI3DWithTraceTable.registered_trace_session1(i,1); %get the registered struct of current ROI

    if ~isempty(tempROI.selected_Z)             %if the ROI have registered trace label
        selectedIdx =tempROI.selected_Z;        %assign the selected registered trace label
        tempROITrace = tempROI.(selectedIdx);   
        ROICorrIdx(i,2) = tempROITrace;         %store the trace of the ROI
    else
        ROICorrIdx(i,2) = {0};                  %if the ROI have no registered trace label, store 0
    end
end

% Remove the rows where the second cell element is 0
ROICorrIdx(cellfun(@(x) isequal(x, 0), ROICorrIdx(:,2)), :) = []; 

pairROINum = size(ROICorrIdx,1);        %get the number of ROIs with value
ROIPairMat = zeros(pairROINum); %create a cell array to store the index of ROI and the corresponding trace

%% calculate the correlation coefficient of ROIs and show the plot
wb = waitbar(0, 'Calculating the correlation coefficient of ROIs...'); %create a waitbar to show the progress
for i = 1 : pairROINum^2
    waitbar(i/(pairROINum^2), wb); %update the waitbar
    [row, col] = ind2sub([pairROINum, pairROINum], i); %get the row and column index of the pair
    if row ~= col
        ROIPairMat(row, col) = corr(cell2mat(ROICorrIdx(row,2)), cell2mat(ROICorrIdx(col,2))); %calculate the correlation coefficient
    elseif row == col
        ROIPairMat(row, col) = 1; %set the diagonal element to 1
    end
end
close(wb); %close the waitbar
figure; 
imagesc(ROIPairMat); %plot the correlation coefficient matrix
colorbar; %add the colorbar
title('Correlation Coefficient Matrix of ROIs'); %add the title
