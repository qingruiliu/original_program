%% step1 load the manually adjusted 3D ROI .tif file from cellPose
%the feature-added version of this program is used.
[fileName,path] = uigetfile('.tif');
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

%%
fillWB = waitbar(0,'0/0','Name','filling the empty ROIs...');
filledVolume = zeros(size(cellPoseVolume));

for label = 1 : length(roiLabels)
    wbStr = append(num2str(label), ' / ' ,num2str(length(roiLabels)));
    waitbar(label/length(roiLabels),fillWB,wbStr);
    %binary mask for current ROI outline
    roiMask = double(cellPoseVolume == roiLabels(label));

    

    for sliceN = 1 : size(filledVolume,3)
        filledsliceN = imfill(roiMask(:,:,sliceN),'holes');
        filledVolume(:,:,sliceN) = filledVolume(:,:,sliceN) + filledsliceN * roiLabels(label);
    end
    currentROIFP = double(filledVolume == label);
    ROI3DWithTraceTable.FP_3D(label) = {currentROIFP};  %save the 3D mask as a 3D matrix of each ROI
end
close(fillWB)

preview = vol3d('CData',filledVolume,'texture','3D');
view(3);

%% 将每一层保存为用于cellReg的输入，并将每个其中每个3D ROI的序号填入table中

%create coordinate table
tableTitle = {'index','planeN','ROINum','ROIIdxList','ROIfootprints','cellRegInput'};
tableTitleTypes = {'double','double','double','cell','cell','cell'};
coordinateTable = table('size',[sliceN length(tableTitle)],...
                  'VariableTypes',tableTitleTypes,'VariableNames',tableTitle);
coordinateTable.index(:) = 1:sliceN;
clearvars tableTitle tableTitleTypes  

[~,name,~] = fileparts(fileName);
startZ = str2double(extractBefore(name,'-'));  %start Z plane number
endZ = startZ + sliceN;          %end Z plane number


for i = 1 : sliceN
    %the 1st z-plane is truncated since python start from 0 and MATLAB is 1
    tempPlane = filledVolume(:,:,i);
    tempROIList = unique(tempPlane);
    tempROIList(tempROIList == 0) = [];            %exclude zero
    
    %write the values to the table elements
    coordinateTable.planeN(i) = startZ+i;  
    coordinateTable.ROINum(i) = nnz(tempROIList);  %exclude zero
    coordinateTable.ROIIdxList(i) = {tempROIList'};

    tempCellRegInput = zeros(512,512,length(tempROIList));
    tempCell = cell(length(tempROIList),1);           %save the temp spatial footprint
    
    for j = 1 : length(tempROIList)
        ROILinearFP = find(tempPlane == tempROIList(j));  %the ROI equals to the j-TH ROI idx
        [x,y] = ind2sub([512 512], ROILinearFP);
        tempCell{j} = [x,y];
        tempCellRegInput(:,:,j) = double(tempPlane == tempROIList(j));
    end

    coordinateTable.ROIfootprints(i) = {tempCell};
    coordinateTable.cellRegInput(i) = {permute(tempCellRegInput,[3 2 1])};
end

save('coordinateTable.mat','coordinateTable','-v7.3')
save('ROI3DWithTraceTable.mat','ROI3DWithTraceTable','-v7.3')

%% check the result (no problem)

Z50 = filledVolume(:,:,50);
FPZ50 = coordinateTable.ROIfootprints{50};
blankImg = zeros(512,512);
for i = 1 : length(FPZ50)
    idxMat = FPZ50{i};
    for j = 1: length(idxMat)
    blankImg(idxMat(j,1),idxMat(j,2)) = 1;
    end
end

%% save the individual cellReg input as individual .mat file

for i = 1 : coordinateTable.index(end)

    cellRegInput = coordinateTable.cellRegInput{i};
    tempSaveStr = append('Z',num2str(coordinateTable.planeN(i)),'.mat');
    save(tempSaveStr,"cellRegInput");
end
