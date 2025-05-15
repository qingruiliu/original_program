%% step1 load the manually adjusted 3D ROI .tif file from cellPose
%the feature-added version of this program is used.
questdlg('Please select the 3D ROI .tif file from cellPose','3D ROI file selection','OK','OK');
[fileName,path] = uigetfile('*masks.tif');
cd(path)
cellPoseVolume = double(tiffreadVolume(fileName));

%% step2(optional):如果cellpose中的mask不是填充好的情况下，需要将每一个区域进行填充，
%  并新建一个3D ROI footprint table

%get the number list of all of the ROIs 
roiLabels = unique(cellPoseVolume);
roiLabels(roiLabels == 0) = [];

%create a blank table of 3D ROIs

tableTitle = {'ROI3DIdx','FP_3D','S1_registered','S2_registered','S3_registered'};
tableTitleTypes = {'double','cell','struct','struct','struct'};
ROI3DWithTraceTable = table('size',[length(roiLabels) length(tableTitle)],...
                            'VariableTypes',tableTitleTypes,'VariableNames',tableTitle);
ROI3DWithTraceTable.ROI3DIdx(:) = 1 : length(roiLabels);
clearvars tableTitle tableTitleTypes

%% extract the 3D ROIs and save the filled ROIs to the table
fillWB = waitbar(0, '0/0', 'Name', 'Filling the empty ROIs...');
filledVolume = zeros(size(cellPoseVolume)); % 初始化填充后的 2D 矩阵

sliceN = size(cellPoseVolume, 3); % 获取切片数量

% 遍历每一层进行处理
for sliceIdx = 1:sliceN
    wbstr = append('Slice ', num2str(sliceIdx), '/', num2str(sliceN));
    waitbar(sliceIdx / sliceN, fillWB, wbstr);

    % 获取当前切片
    currentSlice = cellPoseVolume(:, :, sliceIdx);

    % 获取当前切片中的 ROI 标签
    roiLabels = unique(currentSlice);
    roiLabels(roiLabels == 0) = []; % 排除背景

    % 初始化当前切片的填充结果
    filledSlice = zeros(size(currentSlice));

    % 遍历当前切片中的每个 ROI
    for labelIdx = 1:length(roiLabels)
        roiMask = (currentSlice == roiLabels(labelIdx)); % 获取当前 ROI 掩码
        filledROI = imfill(roiMask, 'holes');           % 填充当前 ROI
        filledSlice(filledROI) = roiLabels(labelIdx);   % 保存填充结果
    end

    % 将填充后的切片保存到 filledVolume
    filledVolume(:, :, sliceIdx) = filledSlice;
end

close(fillWB);

preview = vol3d('CData',filledVolume,'texture','3D'); %preview the filled volume
view(3);
%delete the unnecessary variables 
clearvars -except filledVolume cellPoseVolume fileName path roiLabels

%% filepath: /home/liu/github/original_program/3D_registration_pipeline/step4_cellPoseSliceAnd3DROI.m

% 创建 coordinateTable
sliceN = size(filledVolume, 3);
tableTitle = {'index', 'planeN', 'ROINum', 'ROIIdxList', 'ROIfootprints', 'cellRegInput'};
tableTitleTypes = {'double', 'double', 'double', 'cell', 'cell', 'cell'};
coordinateTable = table('size', [sliceN, length(tableTitle)], ...
                        'VariableTypes', tableTitleTypes, 'VariableNames', tableTitle);
coordinateTable.index(:) = 1:sliceN;
clearvars tableTitle tableTitleTypes  

%% 提取切片编号
startZ = 0;  % 提取文件名中的数字
endZ = 1+sliceN;  % 计算结束的 Z 平面编号

% 初始化进度条
wb1 = waitbar(0, '0/0', 'Name', 'Saving the 3D ROI data...');

% 分块保存 coordinateTable 的数据
for i = 1:sliceN
    waitbar(i / sliceN, wb1, append(num2str(i), '/', num2str(sliceN)));
    tempPlane = filledVolume(:, :, i);
    tempROIList = unique(tempPlane);
    tempROIList(tempROIList == 0) = [];  % 排除背景

    % 写入表格
    coordinateTable.planeN(i) = i;
    coordinateTable.ROINum(i) = nnz(tempROIList);
    coordinateTable.ROIIdxList(i) = {tempROIList'};

    % 如果没有 ROI，赋值为空
    if isempty(tempROIList)
        coordinateTable.ROIfootprints(i) = {{}};
        coordinateTable.cellRegInput(i) = {[]};
        continue;
    end

    tempCellRegInput = zeros(size(filledVolume, 1), size(filledVolume, 2), length(tempROIList));
    tempCell = cell(length(tempROIList), 1);

    for j = 1:length(tempROIList)
        roiMask = (tempPlane == tempROIList(j));
        [x, y] = find(roiMask);  % 获取 ROI 的坐标
        tempCell{j} = [x, y];
        tempCellRegInput(:, :, j) = double(roiMask);
    end

    coordinateTable.ROIfootprints(i) = {tempCell};
    coordinateTable.cellRegInput(i) = {permute(tempCellRegInput, [3, 2, 1])};

end
close(wb1);

% 保存 ROI3DWithTraceTable
save('ROI3DWithTraceTable.mat', 'ROI3DWithTraceTable', '-v7.3');

%% 保存 individual cellReg input
questdlg('Please select the folder to save the sliced cellPose ROI', 'cellReg input saving', 'OK', 'OK');
path = uigetdir();
cd(path);

for i = 1:coordinateTable.index(end)
    cellRegInput = coordinateTable.cellRegInput{i};
    tempSaveStr = append('Z', num2str(coordinateTable.planeN(i)), '.mat');
    save(tempSaveStr, "cellRegInput");
end

clearvars -except *Table
%% check the result (no problem)

%Z50 = filledVolume(:,:,50);           
%FPZ50 = coordinateTable.ROIfootprints{50}; % 
%blankImg = zeros(512,512);
%for i = 1 : length(FPZ50)
%    idxMat = FPZ50{i};
%    for j = 1: length(idxMat)
%    blankImg(idxMat(j,1),idxMat(j,2)) = 1;
%    end
%end

%% save the individual cellReg input as individual .mat file
questdlg('Please select the folder to save the sliced cellPose ROI','cellReg input saving','OK','OK');
path = uigetdir();
cd(path)

for i = 1 : coordinateTable.index(end)
    cellRegInput = coordinateTable.cellRegInput{i};
    tempSaveStr = append('Z',num2str(coordinateTable.planeN(i)),'.mat');
    save(tempSaveStr,"cellRegInput");
end

clearvars -except *Table