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

% Remove the ROIs without trace
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
        %only calculate the distance of x
        ROIPairDistantMat(row, col) = norm(ROICorrIdx{row,3}(:,1:2) - ROICorrIdx{col,3}(:,1:2));
    elseif row == col  %the same ROIs
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
xlim([0 100]);
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

%read the volume of the tiff file
questdlg('Please select the volume imaging of current mouse','Select the tiff file of the volume','OK','OK');
[volumeName, volumePath] = uigetfile('*.tif','Select the tiff file of the volume');
cd(volumePath);
dendriteVolume = double(tiffreadVolume(volumeName))/255;
%volumeViewer(dendriteVolume); %show the volume of the L5 neuron
%% 
questdlg('proceed to the next step','Next Step','OK','OK');

denoiseNet = denoisingNetwork('DnCNN'); %create a denoising network
apicalDendriteMat = zeros(size(dendriteVolume,1),size(dendriteVolume,2),50); %create a matrix to store the apical dendrite orientation
orientationTable = zeros(5,3); %the table to save the orientation information of each 10 slices

wb3 = waitbar(0, 'Processing the volume...'); %create a waitbar to show the progress
for i = 1:50    %depth 150 - 300
    waitbar(i/50, wb3); %update the waitbar
    tempZ = dendriteVolume(:,:,i);
    tempZdenoise1 = denoiseImage(tempZ,denoiseNet); %denoise the image

    tempZdenoise = imbinarize(tempZdenoise1,0.5); %binarize the image
    
    %remove the large connected component
    tempZdenoise = bwareafilt(tempZdenoise,[2 30]); %remove the large connected component
    apicalDendriteMat(:,:,i) = tempZdenoise;
    if mod(i,10) == 0
        tempVolume = apicalDendriteMat(:,:,i-9:i);
        orientationTable(i/10,:) = table2array(regionprops3(tempVolume,'Orientation'));
    end
end
close(wb3); %close the waitbar
volumeViewer(apicalDendriteMat); %show the volume of the apical dendrite

%% check the connected component of the apical dendrite
connectedApicalDendrite = bwconncomp(apicalDendriteMat,18); %get the connected component of the apical dendrite
%exclude the small connected component
for i = 1 : connectedApicalDendrite.NumObjects
    if length(connectedApicalDendrite.PixelIdxList{i}) < 10   %remove the dendrite elements smaller than 10 pixels
        connectedApicalDendrite.PixelIdxList{i} = {};
    end
end
%remove the empty cell
connectedApicalDendrite.PixelIdxList = connectedApicalDendrite.PixelIdxList(~cellfun('isempty',connectedApicalDendrite.PixelIdxList));
%change the value of NumObjects
connectedApicalDendrite.NumObjects = length(connectedApicalDendrite.PixelIdxList);
%% since the z-step of each slice is 3um, the the elongated matrix with 3um step
% 'nearest or linear'
apicalDendriteMatwithZstep = imresize3(apicalDendriteMat,[size(apicalDendriteMat,1),size(apicalDendriteMat,2),size(apicalDendriteMat,3)*3],'linear');
volumeViewer(apicalDendriteMatwithZstep); %show the volume of the apical dendrite with 3um step

 
%% calculate the orientation of each apical dendrite
totalOrientation = regionprops3(connectedApicalDendrite,'EigenVectors','EigenValues'); %get the orientation of the apical dendrite

%calculate the mean eigenvectors and eigenvalues of the apical dendrite
meanEigenVectors = zeros(3,3);   %create a matrix to store the mean eigenvectors
meanEigenValues = zeros(3,1);    %create a vector to store the mean eigenvalues
totalWeight = 0;                 %initialize total weight
for i = 1 : connectedApicalDendrite.NumObjects
    weight = norm(totalOrientation.EigenValues{i}); %use the norm of eigenvalues as weight
    meanEigenVectors = abs(meanEigenVectors) + weight * totalOrientation.EigenVectors{i};
    meanEigenValues = meanEigenValues + weight * totalOrientation.EigenValues{i};
    totalWeight = totalWeight + weight;             %accumulate the total weight
end
meanEigenVectors = meanEigenVectors / totalWeight;  %weighted mean eigenvectors
meanEigenValues = meanEigenValues / totalWeight;    %weighted mean eigenvalues

%% transform the 3D ROI center using the mean eigenvectors and eigenvalues
%ROICorrIdx(:,4) = (); %create a cell array to store the transformed 3D ROI center
for i = 1 : size(ROICorrIdx,1)
    tempCenter = cell2mat(ROICorrIdx(i,3));
    normEigenValueMat = diag(meanEigenValues / norm(meanEigenValues));
    transformedCenter = (meanEigenVectors * normEigenValueMat * tempCenter')';   %transform the 3D ROI center
    ROICorrIdx(i,4) = {transformedCenter}; %store the transformed 3D ROI center
end

%% code from other files
%NVec: unit vector representing the microcolumn axis or apical dendrites
ClmAxisVec=NVec; 

ClmAxisVec(3)=sqrt(1-(ClmAxisVec(1)^2+ClmAxisVec(2)^2));
[tgtrng,tgtang] = rangeangle(ClmAxisVec');

TempAngle=[tgtang(1),tgtang(2)-90];

Theta2=deg2rad(-1*TempAngle(1));
Phi2=deg2rad(-1*TempAngle(2));
ClmPrjMtx1=[cos(Theta2) -sin(Theta2);sin(Theta2) cos(Theta2)];
ClmPrjMtx2=[cos(Phi2) -sin(Phi2);sin(Phi2) cos(Phi2)];

TempVec=ClmAxisVec;

Temp4=ClmPrjMtx1*[TempVec(1);TempVec(2)];
TempVec(1)=Temp4(1);
TempVec(2)=Temp4(2);
Temp3=ClmPrjMtx2*[TempVec(1);TempVec(3)];
TempVec(1)=Temp3(1);
TempVec(3)=Temp3(2);
[temprng,tempang] = rangeangle(TempVec');

disp('ClmAxisVecの回転後の行列')
TempVec

%%%

TempPosUM1=Pos; %Pos: 3D coordinates of GCaMP6 cells


tTempPosUM1=TempPosUM1;
h = waitbar(0,'--座標回転中1--');
for iCell=1:size(TempPosUM1,1)
    Temp1=ClmPrjMtx1*TempPosUM1(iCell,1:2)';
    tTempPosUM1(iCell,1)=Temp1(1);
    tTempPosUM1(iCell,2)=Temp1(2);
    Temp2=ClmPrjMtx2*[tTempPosUM1(iCell,1);TempPosUM1(iCell,3)];
    tTempPosUM1(iCell,1)=Temp2(1);
    tTempPosUM1(iCell,3)=Temp2(2);

    waitbar(iCell/size(tTempPosUM1,1));
end
close(h)