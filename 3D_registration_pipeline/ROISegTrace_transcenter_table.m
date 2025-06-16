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
        ROISegIdx(i,2) = {tempROITrace};         %store the trace of the ROI
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
ROISegTraceTable = cell2table(ROISegIdx,'VariableNames',{'ROIIndex','Segmented_trace','Center','TransformedCenter'});

%% use the PCA method to plot the fitting plane as the tangential plane
centroidsMat = cell2mat(ROISegIdx(:,3));
centroidsX = centroidsMat(:,1);
centroidsY = centroidsMat(:,2);
centroidsZ = centroidsMat(:,3);

% Perform PCA on the centroid data to find the major component
[coeff,~,~] = pca(centroidsMat);

normalVector = coeff(:,3); %get the normal vector of the fitting plane
meanPoint = mean(centroidsMat,1); %get the mean point of the fitting plane

%visualize the data and the fitting plane
figure;
scatter3(centroidsX, centroidsY, centroidsZ,30,[0.3 0.3 0.3],'filled');
hold on;

%define the plane given the PCA normal vector
[xPlane, yPlane] = meshgrid(linspace(min(centroidsX),max(centroidsX),10),linspace(min(centroidsY),max(centroidsY),10));

zPlane = meanPoint(3) - normalVector(1)/normalVector(3)*(xPlane - meanPoint(1)) - normalVector(2)/normalVector(3)*(yPlane - meanPoint(2)); %calculate the Z value of the plane

%plot the fitting plane
surf(xPlane, yPlane, zPlane,'FaceAlpha',0.3,'EdgeColor','none');
colormap('abyss');

xlabel('X');ylabel('Y');zlabel('Z');
%title('Fitting Plane of the ROIs');
grid on;
axis equal;
xlabel('Original X (μm)')
xticks(0:100:500)
ylabel('Original Y (μm)')
yticks(0:100:500)
zlabel('Depth (μm)')
zticks(10:50:140)
zticklabels({'360','400','440'})
set(gca,'ZDir','reverse')
set(gca,'Fontsize',16)
%legend('Data Points','Fitting Plane');
view(3);
grid off
set(gca,'FontSize',20)
hold off;

%display the normal vector and the angles with the x,y,z axis
disp('Normal Vector of the Fitting Plane:');
disp(normalVector);

%calculate the angles between the normal vector and the x,y,z axis
angleXY = acosd(abs(normalVector(3)/norm(normalVector)));
angleXZ = acosd(abs(normalVector(2)/norm(normalVector)));
angleYZ = acosd(abs(normalVector(1)/norm(normalVector)));

%transformed the centroids based on PCA result
centeredCentroids = centroidsMat - meanPoint;
transMat = coeff;
transformedCentroids = centeredCentroids * transMat;

XYplaneData = transformedCentroids(:,1:2); %get the data of the XY plane

%visualize the points in the new XY plane
figure;
scatter(XYplaneData(:,1), XYplaneData(:,2),'k','filled');
xlabel('X');ylabel('Y');
title('Centroids in the transformed XY Plane');
grid on; 
axis equal;

ROISegTraceTable.TransformedCenter = transformedCentroids;
savestr = append('ROISegTraceTable_',sessionNum);
save(savestr,'ROISegTraceTable','-v7.3')