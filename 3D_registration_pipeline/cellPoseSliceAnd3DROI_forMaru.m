%% load the manually adjusted 3D ROI .tif file from cellPose
%the feature-added version of this program is used.
questdlg('Please select the 3D ROI .tif file from cellPose','3D ROI file selection','OK','OK');
[fileName,path] = uigetfile('*masks.tif');
cd(path)
cellPoseVolume = double(tiffreadVolume(fileName));

%%  create a new table
%get the number list of all of the ROIs 
roiLabels = unique(cellPoseVolume);
roiLabels(roiLabels == 0) = [];

%create a blank table of 3D ROIs
tableTitle = {'ROI3DIdx','FP_3D','centroid'};
tableTitleTypes = {'double','cell','cell'};
ROI3DTable = table('size',[length(roiLabels) length(tableTitle)],...
                            'VariableTypes',tableTitleTypes,'VariableNames',tableTitle);
ROI3DTable.ROI3DIdx(:) = 1 : length(roiLabels);
clearvars tableTitle tableTitleTypes

%% extract the 3D ROIs and save the filled ROIs to the table
fillWB = waitbar(0, '0/0', 'Name', 'Filling the empty ROIs...');
filledVolume = zeros(size(cellPoseVolume)); %blank filled volume

sliceN = size(cellPoseVolume, 3); %z-depth

for sliceIdx = 1:sliceN
    wbstr = append('Slice ', num2str(sliceIdx), '/', num2str(sliceN));
    waitbar(sliceIdx / sliceN, fillWB, wbstr);

    % current slice
    currentSlice = cellPoseVolume(:, :, sliceIdx);

    % ROI labels on current slice
    roiLabel2 = unique(currentSlice);
    roiLabel2(roiLabel2 == 0) = [];

    filledSlice = zeros(size(currentSlice));

    for labelIdx = 1:length(roiLabel2)
        roiMask = (currentSlice == roiLabel2(labelIdx)); % current ROI label
        filledROI = imfill(roiMask, 'holes');           % fill current ROI
        filledSlice(filledROI) = roiLabel2(labelIdx);   % save into filled slice
    end

    % fill the current slice in filledVolume
    filledVolume(:, :, sliceIdx) = filledSlice;
end

close(fillWB);
%save the filledVolume
save('filledVolume.mat','filledVolume')
preview = vol3d('CData',filledVolume,'texture','3D'); %preview the filled volume
view(3);

%% build the 3D ROI table
wb = waitbar(0, '0/0', 'Name', 'Building the 3D ROI table...');
for i = 1:length(roiLabels)
    wbstr = append(num2str(i), '/', num2str(length(roiLabels)));
    waitbar(i / length(roiLabels), wb, wbstr);
    
    % Get the indices of the current ROI in the 3D volume
    [x, y, z] = ind2sub(size(filledVolume), find(filledVolume == roiLabels(i)));
    
    % Save the coordinates as a cell array of [x, y, z]
    ROI3DTable.FP_3D{i} = [x, y, z];
    ROI3DTable.centroid{i} = mean([x,y,z]);
end
save('ROI3DTable.mat','ROI3DTable','-v7.3')
close(wb)

%% Function to visualize a selected ROI in the original 3D volume
selectedROI = inputdlg('Enter the ROI index to visualize:', 'Select ROI', [1 50]);
selectedROI = str2double(selectedROI);

if ~isnan(selectedROI) && ismember(selectedROI, roiLabels)
    % Get the coordinates of the selected ROI
    roiCoords = ROI3DTable.FP_3D{selectedROI};
    
    % Create a blank volume to display only the selected ROI
    roiVolume = zeros(size(cellPoseVolume));
    for k = 1:size(roiCoords, 1)
        roiVolume(roiCoords(k, 1), roiCoords(k, 2), roiCoords(k, 3)) = cellPoseVolume(roiCoords(k, 1), roiCoords(k, 2), roiCoords(k, 3));
    end
    
    % Preview the volume containing only the selected ROI
    preview = vol3d('CData', roiVolume, 'texture', '3D');
    view(3);
    title(['Visualization of ROI ', num2str(selectedROI)]);
else
    msgbox('Invalid ROI index selected.', 'Error', 'error');
end

%% 创建 coordinateTable
sliceN = size(filledVolume, 3);
tableTitle = {'index', 'planeN', 'ROINum', 'ROIIdxList', 'ROIfootprints'};
tableTitleTypes = {'double', 'double', 'double', 'cell', 'cell'};
coordinateTable = table('size', [sliceN, length(tableTitle)], ...
                        'VariableTypes', tableTitleTypes, 'VariableNames', tableTitle);
coordinateTable.index(:) = 1:sliceN;
clearvars tableTitle tableTitleTypes  

% get slice number
wb1 = waitbar(0, '0/0', 'Name', 'Saving the 3D ROI data...');

for i = 1:sliceN
    waitbar(i / sliceN, wb1, append(num2str(i), '/', num2str(sliceN)));
    tempPlane = filledVolume(:, :, i);
    tempROIList = unique(tempPlane);
    tempROIList(tempROIList == 0) = []; 

    % write in the table
    coordinateTable.planeN(i) = i;
    coordinateTable.ROINum(i) = nnz(tempROIList);
    coordinateTable.ROIIdxList(i) = {tempROIList'};

    % if no ROIs in current slice
    if isempty(tempROIList)
        coordinateTable.ROIfootprints(i) = {{}};
        continue;
    end

    tempCell = cell(length(tempROIList), 1);

    for j = 1:length(tempROIList)
        roiMask = (tempPlane == tempROIList(j));
        [x, y] = find(roiMask);  % get the coordinates
        tempCell{j} = [x, y];
        %tempCellRegInput(:, :, j) = double(roiMask);
    end

    coordinateTable.ROIfootprints(i) = {tempCell};

end
close(wb1);
save('coordinateTable.mat', 'coordinateTable', '-v7.3');


