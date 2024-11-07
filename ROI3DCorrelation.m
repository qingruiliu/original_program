%% index back to the whole trace 
ROICorrIdx = cell(height(ROI3DWithTraceTable),3);   %create a cell array to store the index of ROI and the corresponding trace

wb1 = waitbar(0, 'Indexing the ROI and trace...'); %create a waitbar to show the progress
for i = 1 : height(ROI3DWithTraceTable)
    waitbar(i/height(ROI3DWithTraceTable), wb1); %update the waitbar
    ROICorrIdx(i,1) = {i};                          %store the index of ROI
    tempROI = ROI3DWithTraceTable.registered_trace_session1(i,1); %get the registered struct of current ROI

    if ~isempty(tempROI.selected_Z)             %if the ROI have registered trace label
        selectedIdx =tempROI.selected_Z;        %assign the selected registered trace label
        tempROITrace = tempROI.(selectedIdx);   
        ROICorrIdx(i,2) = tempROITrace;         %store the trace of the ROI
    else
        ROICorrIdx(i,2) = {0};                  %if the ROI have no registered trace label, store 0
    end

    %assign and calculate the center of each ROI
    tempFP_3D = ROI3DWithTraceTable.FP_3D{i};   %get the footprint of current ROI
    [x,y,z] = size(tempFP_3D);                   %get the size of the footprint
    [roi_x, roi_y,roi_z] = ind2sub([x,y,z],find(tempFP_3D == 1));     %get the index of the footprint
    tempCenter = [mean(roi_x), mean(roi_y), mean(roi_z)];              %calculate the center of the footprint
    ROICorrIdx(i,3) = {tempCenter};                                     %store the center of the ROI
end
close(wb1); %close the waitbar

% Remove the rows where the second cell element is 0
ROICorrIdx(cellfun(@(x) isequal(x, 0), ROICorrIdx(:,2)), :) = []; 

pairROINum = size(ROICorrIdx,1);        %get the number of ROIs with value
ROIPairCorrMat = zeros(pairROINum); %create a cell array to store the index of ROI and the corresponding trace
ROIPairDistantMat = zeros(pairROINum); %create a cell array to store the index of ROI and the corresponding trace
%% calculate the correlation coefficient of ROIs and show the plot

wb2 = waitbar(0, 'Calculating the correlation coefficient of ROIs...'); %create a waitbar to show the progress
%plot the cell-pair distance and correlation coefficient on the same figure
for i = 1 : pairROINum^2
    waitbar(i/(pairROINum^2), wb2); %update the waitbar
    [row, col] = ind2sub([pairROINum, pairROINum], i); %get the row and column index of the pair
    if row ~= col      %different ROIs
        ROIPairCorrMat(row, col)  = corr(cell2mat(ROICorrIdx(row,2)), cell2mat(ROICorrIdx(col,2))); %calculate the correlation coefficient
        ROIPairDistantMat(row, col) = norm(cell2mat(ROICorrIdx(row,3)) - cell2mat(ROICorrIdx(col,3))); %calculate the distance between the centers of the ROIs 
    elseif row == col  %same ROI
        ROIPairCorrMat(row, col) = 1; %set the diagonal element to 1
        ROIPairDistantMat(row, col) = 0; %set the diagonal element to 0
    end
end
close(wb2); %close the waitbar
%% plot the correlation and distance matrix
figure; 
imagesc(ROIPairCorrMat); %plot the correlation coefficient matrix
colorbar; %add the colorbar
clim([-0.2 0.2]);
title('Correlation Coefficient Matrix of ROIs'); %add the title
set(gca,'FontSize',16)
figure;
imagesc(ROIPairDistantMat); %plot the distance matrix
colorbar; %add the colorbar
clim([0 100]);
title('Distance Matrix of ROIs'); %add the title
set(gca,'FontSize',16)

%% plot the correlation coefficient and distance of ROIs on the same figure
figure;
hold on
xlim([0 inf]);
ylim([-0.2 1]);

corrAndDist(:,1) = ROIPairCorrMat(:);   %reshape the correlation coefficient matrix to a column vector
corrAndDist(:,2) = ROIPairDistantMat(:);%reshape the distance matrix to a column vector

%sort the correlation coefficient and distance matrix by distance
corrAndDist = sortrows(corrAndDist,2);
%exclude the diagonal element
corrAndDist1 = corrAndDist(corrAndDist(:,2) ~= 0,:);
scatter(corrAndDist1(:,2),corrAndDist1(:,1),10,'Filled'); %plot the scatter plot of correlation coefficient and distance
title('Correlation Coefficient and Distance of ROIs'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndDist1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndDist(:,2) >= distanceBin(i) & corrAndDist(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndDist(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(1:end-1)+5,corrBin,'k','LineWidth',2);  %plot the average correlation coefficient of the distance bin
hold off
xlabel('Distance between the centers of ROI pairs');     %add the x-axis label
ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
set(gca,'FontSize',16)

%% adjust the direction to the actual tangential and radial orientation