%% re-arrange the output of suite2p into MATLAB table format

% select the target Fall.mat file
disp('Select the ***Fall.mat*** of current processing imaging plane')
[fileName,path] = uigetfile('.mat');
cd(path)
load(fileName);
fprintf('%s LOADED \n',path)

%% create a empty table to save different variables

tableTitle = {'ROIindexS2P','isCell','Prob','rawF','DeconvF','Stat','cellRegInputS2P'};
tableTitleTypes = {'double','double','double','cell','cell','cell','cell'};

%create the empty table
suite2pTable = table('size',[length(iscell) length(tableTitle)],...
                        'VariableNames',tableTitle, ...
                        'VariableTypes',tableTitleTypes);

%clearvars tableTitleTypes tableTitle

%% build the suite2pTable

suite2pTable.isCell(:) = iscell(:,1);
suite2pTable.Prob(:) = iscell(:,2);

for i = 1 : length(iscell)
    emptyImg = zeros(512,512);
    suite2pTable.ROIindexS2P(i) = i;
    suite2pTable.rawF(i) = {F(i,:)'};
    suite2pTable.DeconvF(i)= {spks(i,:)'};
    suite2pTable.Stat(i) = stat(i);
    tempX = double(stat{i}.xpix)';
    tempY = double(stat{i}.ypix)';
    tempLam = stat{i}.lam';

    %max-min normalize the pixel contribution of Lam
    tempLamNorm = (tempLam - min(tempLam)) / (max(tempLam) - min(tempLam));
    tempSoma = stat{i}.soma_crop';
    emptyImg(sub2ind([512 512],tempY(tempSoma),tempX(tempSoma))) = tempLamNorm(tempSoma); 
    suite2pTable.cellRegInputS2P(i) = {emptyImg};
end

%delete all the trace if not cell
% toDelete = suite2pTable.isCell == 0;
% suite2pTable(toDelete,:) = [];
% suite2pTable.ROIindex = (1:size(suite2pTable,1))';

% sort the rows based on the isCell and ROIindex
suite2pTable = sortrows(suite2pTable,{'isCell','ROIindexS2P'},{'descend','ascend'});
save('suite2pTable.mat','suite2pTable')
clearvars -except *Table path

% plot the ROI again and compare with the original image (或者手动加载suite2pTable.mat之后)

suite2pImage = zeros(512,512);
for i = 1 : size(suite2pTable)
    if suite2pTable.isCell(i) == 1  %only get the isCell spatial footprint
    suite2pImage = suite2pImage + suite2pTable.cellRegInputS2P{i,1};
    end
end
figure; 
imshow(suite2pImage)

title('spatial footprint of suite2p, isCell only','FontSize',20)

% get the iscell 3D matrix for cellReg input, save to cellRegInput.mat

isCellNum = sum(suite2pTable.isCell);
suite2pInput = zeros(512,512,isCellNum);

for i = 1 : isCellNum
    suite2pInput(:,:,i) = suite2pTable.cellRegInputS2P{i};
end

%% 
cellRegInput = permute(suite2pInput,[3 2 1]);
saveStr = append('cellRegInputZ',extractAfter(path,'Z'));
saveStr = append(extractBefore(saveStr,'/'),'.mat');
save(saveStr,"cellRegInput");

clearvars -except *Table