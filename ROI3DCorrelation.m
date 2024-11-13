%% adjust the direction to the actual tangential and radial orientation

%read the volume of the tiff file
questdlg('Please select the volume imaging of current mouse','Select the tiff file of the volume','OK','OK');
[volumeName, volumePath] = uigetfile('*.tif','Select the tiff file of the volume');
cd(volumePath);
dendriteVolume = double(tiffreadVolume(volumeName))/255;
%volumeViewer(dendriteVolume); %show the volume of the L5 neuron

%% reconstruct the volume of the apical dendrite
questdlg('proceed to the next step','Next Step','OK','OK');

denoiseNet = denoisingNetwork('DnCNN'); %create a denoising network
apicalDendriteMat = zeros(size(dendriteVolume,1),size(dendriteVolume,2),50); %create a matrix to store the apical dendrite orientation

wb3 = waitbar(0, 'Processing the volume...'); %create a waitbar to show the progress
for i = 1:50    %depth 150 - 300
    waitbar(i/50, wb3); 
    tempZ = dendriteVolume(:,:,i);
    tempZdenoise1 = denoiseImage(tempZ,denoiseNet); %denoise the image

    tempZdenoise = imbinarize(tempZdenoise1); %binarize the image
    
    %remove the large connected component
    tempZdenoise = bwareafilt(tempZdenoise,[0 25]); %remove the large connected component
    apicalDendriteMat(:,:,i) = tempZdenoise;
end
close(wb3); 
apicalDendriteMat = double(bwareaopen(apicalDendriteMat,200,18)); %remove the small connected component less than 100 pixels
%resize the volume with 3-times z-step
apicalDendriteMat = imresize3(apicalDendriteMat,[size(apicalDendriteMat,1),size(apicalDendriteMat,2),size(apicalDendriteMat,3)*3],'linear');

volumeViewer(apicalDendriteMat); %show the volume of the apical dendrite
%% check the connected component of the apical dendrite
connectedApicalDendrite = bwconncomp(apicalDendriteMat,18); 

% calculate the orientation of each apical dendrite
totalProps = regionprops3(connectedApicalDendrite, 'Orientation', 'Centroid'); % get the orientation and principal axis length of the apical dendrite

%plot the centroid and the orientation of each apical dendrite in a 3d volume
figure;
x_quiver = totalProps.Centroid(:,1);
y_quiver = totalProps.Centroid(:,2);
z_quiver = totalProps.Centroid(:,3);

u_quiver = cosd(totalProps.Orientation(:,1)).*sind(totalProps.Orientation(:,2));
v_quiver = sind(totalProps.Orientation(:,1)).*sind(totalProps.Orientation(:,2));
w_quiver = cosd(totalProps.Orientation(:,2));
quiver3(x_quiver, y_quiver, z_quiver, u_quiver, v_quiver, w_quiver,'LineWidth', 1);
axis equal;

% manually adjust the oddball element in the quiver plot

%use the average quiver vector to generate the euler angle
meanVector = [mean(u_quiver), mean(v_quiver), mean(w_quiver)];
% Normalize the mean vector
meanVector = meanVector / norm(meanVector);

% Calculate the rotation matrix using Rodrigues' rotation formula
theta = acos(meanVector(3)); % Angle between meanVector and z-axis
k = [-meanVector(2), meanVector(1), 0]; % Rotation axis (cross product of meanVector and z-axis)
k = k / norm(k); % Normalize the rotation axis

K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0]; % Skew-symmetric matrix of k

rotationMatrix = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
disp('Rotation matrix from mean vector:');
disp(rotationMatrix);


%convert the mean vector to the rotation matrix

%% (alternative) Create rotation matrix from mean orientation
% meanOrientation = mean(totalProps.Orientation,1); %get the mean orientation of the apical dendrite
% theta = deg2rad(meanOrientation(1));
% phi = deg2rad(meanOrientation(2));
% psi = deg2rad(meanOrientation(3));
% 
% Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
% Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
% Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% 
% rotationMatrix = Rz * Ry * Rx;
% disp('rorationMatrix: ');
% disp(rotationMatrix)

%% index back to the whole trace 
ROICorrIdx = cell(height(ROI3DWithTraceTable),4);   %create a cell array to store the index of ROI and the corresponding trace

wb1 = waitbar(0, 'Indexing the ROI and trace...'); %create a waitbar to show the progress
for i = 1 : height(ROI3DWithTraceTable)
    waitbar(i/height(ROI3DWithTraceTable), wb1); %update the waitbar
    ROICorrIdx(i,1) = {i};                       %store the index of ROI
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

    %resize the footprint volume to the same scale as apical dendrite volume (3 times z-step)
    tempFP_3D = imresize3(tempFP_3D,[size(tempFP_3D,1),size(tempFP_3D,2),size(tempFP_3D,3)*3],'linear');

    [x,y,z] = size(tempFP_3D);                   %get the size of the footprint
    [roi_x, roi_y,roi_z] = ind2sub([x,y,z],find(tempFP_3D == 1));     %get the index of the footprint, 1 is the positive pixels for this ROI
    tempCenter = [mean(roi_x), mean(roi_y), mean(roi_z)];              %calculate the center of the footprint
    ROICorrIdx(i,3) = {tempCenter};                                     %store the center of the ROI
    transformedCenter = (rotationMatrix * tempCenter')';   %transform the 3D ROI center using rotation matrix
    ROICorrIdx(i,4) = {transformedCenter}; %store the transformed 3D ROI center
end
close(wb1); %close the waitbar

% Remove the ROIs without trace
ROICorrIdx(cellfun(@(x) isequal(x, 0), ROICorrIdx(:,2)), :) = []; 

%save the ROICorrIdx as a table 
ROICorrIdxTable = cell2table(ROICorrIdx,'VariableNames',{'ROIIndex','Trace','Center','TransformedCenter'});
writetable(ROICorrIdxTable,'ROICorrIdxTable.mat'); %save the table as a .mat file

pairROINum = size(ROICorrIdx,1);        %get the number of ROIs with value
ROIPairCorrMat = zeros(pairROINum); %create a matrix to save the correlation coefficient of paired ROIs
ROIPairXYDistMat = zeros(pairROINum); %create a matrix to save the X-Y distance of paired ROIs
ROIPairXYZDistMat = zeros(pairROINum); %create a matrix to save the X-Y-Z distance of paired ROIs
ROIPairTransDistantMat = zeros(pairROINum); %create a matrix to save the transformed X-Y distance of paired ROIs
%% calculate the correlation coefficient of ROIs and show the plot
wb2 = waitbar(0, 'Calculating the correlation coefficient of ROIs...'); %create a waitbar to show the progress
%plot the cell-pair distance and correlation coefficient on the same figure
for i = 1 : pairROINum^2
    waitbar(i/(pairROINum^2), wb2); %update the waitbar
    [row, col] = ind2sub([pairROINum, pairROINum], i); %get the row and column index of the pair
    if row ~= col      %different ROIs
        ROIPairCorrMat(row, col)  = corr(cell2mat(ROICorrIdx(row,2)), cell2mat(ROICorrIdx(col,2))); %calculate the correlation coefficient
        
        ROIPairXYDistMat(row, col) = norm(ROICorrIdx{row,3}(:,1:2) - ROICorrIdx{col,3}(:,1:2));
        ROIPairXYZDistMat(row, col) = norm(ROICorrIdx{row,3} - ROICorrIdx{col,3}); % calculate the distance of the center of the ROIs
        ROIPairTransDistantMat(row, col) = norm(ROICorrIdx{row,4}(:,1:2) - ROICorrIdx{col,4}(:,1:2)); % calculate the distance of the transformed center 
    elseif row == col  %the same ROIs
        ROIPairCorrMat(row, col) = 1; %set the diagonal element to 1
        ROIPairXYDistMat(row, col) = 0; %set the diagonal element to 0
        ROIPairXYZDistMat(row, col) = 0; %set the diagonal element to 0
        ROIPairTransDistantMat(row, col) = 0; %set the diagonal element to 0
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
imagesc(ROIPairXYDistMat); %plot the distance matrix
colorbar; %add the colorbar
clim([0 100]);
title('X-Y Distance Matrix of ROIs'); %add the title
set(gca,'FontSize',16)

figure;
imagesc(ROIPairXYZDistMat); %plot the distance matrix
colorbar; %add the colorbar
clim([0 100]);
title('X-Y-Z Distance Matrix of ROIs'); %add the title
set(gca,'FontSize',16)

figure;
imagesc(ROIPairTransDistantMat); %plot the distance matrix
colorbar; %add the colorbar
clim([0 100]);
title('Transformed distance of ROI pairs'); %add the title
set(gca,'FontSize',16)

%% plot the correlation coefficient and untransformed distance of ROIs on the same figure
figure;
hold on
xlim([0 50]);
ylim([-0.1 0.2]);

corrAndDist(:,1) = ROIPairCorrMat(:);   %reshape the correlation coefficient matrix to a column vector
corrAndDist(:,2) = ROIPairXYDistMat(:);
corrAndDist(:,3) = ROIPairXYZDistMat(:);

%sort the correlation coefficient and distance matrix by distance
corrAndDist = sortrows(corrAndDist,2);

corrAndDist1 = corrAndDist(corrAndDist(:,2) ~= 0,:); %exclude the diagonal element
%exclude the cell pair with XYZ distance smaller than 20 um
corrAndDist1 = corrAndDist1(corrAndDist1(:,3) > 20,:);

scatter(corrAndDist1(:,2),corrAndDist1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('Correlation Coefficient and X-Y Distance of ROIs'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndDist1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndDist1(:,2) >= distanceBin(i) & corrAndDist1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndDist1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(2:end)/2,corrBin,'k','LineWidth',2);  %plot the average correlation coefficient of the distance bin
hold off
xlabel('Distance between the centers of ROI pairs');     %add the x-axis label
ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
set(gca,'FontSize',16)

%% plot the correlation coefficient and transformed distance of ROIs on the same figure
figure;
hold on
xlim([0 50]);
ylim([-0.05 0.1]);

corrAndTransDist(:,1) = ROIPairCorrMat(:);   %reshape the correlation coefficient matrix to a column vector
corrAndTransDist(:,2) = ROIPairTransDistantMat(:);%reshape the distance matrix to a column vector
corrAndTransDist(:,3) = ROIPairXYZDistMat(:);
%sort the correlation coefficient and distance matrix by distance
corrAndTransDist = sortrows(corrAndTransDist,2);
%exclude the diagonal element
corrAndTransDist1 = corrAndTransDist(corrAndTransDist(:,2) ~= 0,:);
%exclude the cell pair with XYZ distance smaller than 20 um
corrAndTransDist1 = corrAndTransDist1(corrAndTransDist1(:,3) > 20,:);

scatter(corrAndTransDist1(:,2),corrAndTransDist1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('Correlation Coefficient and tangential distance of ROIs'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndTransDist1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndTransDist1(:,2) >= distanceBin(i) & corrAndTransDist1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndTransDist1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(2:end)/2,corrBin,'g','LineWidth',2);  %plot the average correlation coefficient of the distance bin

%% generate the surrogate data from toe ROI CorrIdx to shuffle the correspondence between the trace and transformed center
surrogateTimes = 1000; %set the number of surrogate data
surroMat = zeros(size(corrAndTransDist1,1),3,surrogateTimes); %create a matrix to store the surrogate data            %create a matrix to store the surrogate data
%first column: correlation coefficient second column: transformed distance  pages: times of surrogate data
wb4 = waitbar(0, 'Generating the surrogate data...'); %create a waitbar to show the progress
for i = 1 : surrogateTimes
    waitbar(i/surrogateTimes, wb4); %update the waitbar
    for j = 1 : size(corrAndTransDist1,1)

        surroMat(j,1,i) = corr(ROICorrIdx{randi(size(ROICorrIdx,1)),2},ROICorrIdx{randi(size(ROICorrIdx,1)),2}); %randomly shuffle the trace
        surroMat(j,2,i) = norm(ROICorrIdx{randi(size(ROICorrIdx,1)),4}(1:2) - ROICorrIdx{randi(size(ROICorrIdx,1)),4}(1:2)); %randomly shuffle the transformed center
        surroMat(j,3,i) = norm(ROICorrIdx{randi(size(ROICorrIdx,1)),3} - ROICorrIdx{randi(size(ROICorrIdx,1)),3}); %calculate the randomized x-y-z center
    end
    %sort the tempSurro by the transformed distance
    tempSurro = sortrows(surroMat(:,:,i),2);
    %exclude the cell pair with XYZ distance smaller than 20 um
    tempSurro = tempSurro(tempSurro(:,3) > 20,:);
    
    %plot the average corr value of surrogate data in the 5um distance bin
    for k = 1 : length(distanceBin)-1
        tempIdx = find(tempSurro(:,2) >= distanceBin(k) & tempSurro(:,2) < distanceBin(k+1)); %find the index of the correlation coefficient in the distance bin
        corrBin(k) = mean(tempSurro(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
    end
    plot(distanceBin(2:end)/2,corrBin,'Color',[0.3 0.3 0.3 0.15],'LineWidth', 0.5);  %plot the average correlation coefficient of the distance bin
end
close(wb4); %close the waitbar

%% generate the surrogate data and check the significance of the correlation coefficient at different distance bin
%randomly shuffle the correspondence between the correlation coefficient and x-y distance for 1000 times
surrogateTimes = 1000;
corrAndTransDistSurrogate = zeros(size(corrAndTransDist1,1),2,surrogateTimes); %create a matrix to store the surrogate data
wb4 = waitbar(0, 'Generating and plotting the surrogate data...'); %create a waitbar to show the progress
surrogateCorrBinMat = zeros(length(distanceBin)-1,surrogateTimes); %create a matrix to store the average correlation coefficient of each distance bin
for i = 1 : surrogateTimes
    waitbar(i/surrogateTimes, wb4); %update the waitbar
    tempIndex = randperm(size(corrAndTransDist1,1)); %randomly shuffle the index
    corrAndTransDistSurrogate(:,1,i) = corrAndTransDist1(tempIndex,1); %only shuffle the correlation coefficient
    corrAndTransDistSurrogate(:,2,i) = corrAndTransDist1(:,2);       %keep the transformed distance unchanged
    %calculate the average correlation coefficient of the distance bin
    for j = 1 : length(distanceBin)-1
        tempIdx = find(corrAndTransDistSurrogate(:,2,i) >= distanceBin(j) & corrAndTransDistSurrogate(:,2,i) < distanceBin(j+1)); %find the index of the correlation coefficient in the distance bin
        corrBin(j) = mean(corrAndTransDistSurrogate(tempIdx,1,i));  %calculate the average correlation coefficient of the distance bin
    end
    surrogateCorrBinMat(:,i) = corrBin; %store the average correlation coefficient of each distance bin
    %plot the average correlation coefficient of the distance bin
    plot(distanceBin(2:end)/2,corrBin,'Color',[0.3 0.3 0.3 0.1],'LineWidth', 0.5);  %plot the average correlation coefficient of the distance bin
end
close(wb4); %close the waitbar

%plot the highest 2.5% and lowest 2.5% of all the surrogateCorrBinMat
upperBound = prctile(surrogateCorrBinMat,97.5,2); %calculate the upper bound of the surrogate data
lowerBound = prctile(surrogateCorrBinMat,2.5,2); %calculate the lower bound of the surrogate data

plot(distanceBin(2:end)/2,upperBound,'k--','LineWidth',1);  %plot the upper bound of the surrogate data
plot(distanceBin(2:end)/2,lowerBound,'k--','LineWidth',1);  %plot the lower bound of the surrogate data

%% plot the averaged surrogate data and upper/lower 2.5% of the surrogate data
corrAndTransDistSurrogateMean = mean(corrAndTransDistSurrogate,3); %calculate the mean of the surrogate data
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndTransDistSurrogateMean(:,2) >= distanceBin(i) & corrAndTransDistSurrogateMean(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndTransDistSurrogateMean(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(2:end)/2,corrBin,'k','LineWidth',1.5);  %plot the average correlation coefficient of the distance bin

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