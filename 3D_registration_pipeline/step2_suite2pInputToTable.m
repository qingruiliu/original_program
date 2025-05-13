%% re-arrange the output of suite2p into MATLAB table format for multiple imaging planes

% select the main directory containing subfolders for each imaging plane
disp('Select the main directory containing subfolders for each imaging plane')
mainDir = uigetdir();
subDirs = dir(mainDir);
subDirs = subDirs([subDirs.isdir] & ~startsWith({subDirs.name}, '.'));

wb = waitbar(0,'0/0','Name','Processing suite2p data...');
for k = 1:length(subDirs)
    waitbar(k/length(subDirs),wb,append('Processing suite2p data...',num2str(k),'/',num2str(length(subDirs))));
    planePath = fullfile(mainDir, subDirs(k).name);
    disp(['Processing folder: ', planePath]);
    
    % select the target Fall.mat file
    fileName = fullfile(planePath, 'Fall.mat');
    if ~isfile(fileName)
        disp(['Fall.mat not found in ', planePath]);
        continue;
    end
    load(fileName);
    fprintf('%s LOADED \n', planePath);

    %% create a empty table to save different variables

    tableTitle = {'ROIindexS2P','isCell','Prob','rawF','DeconvF','Stat','cellRegInputS2P'};
    tableTitleTypes = {'double','double','double','cell','cell','cell','cell'};

    %create the empty table
    suite2pTable = table('size',[length(iscell) length(tableTitle)],...
                            'VariableNames',tableTitle, ...
                            'VariableTypes',tableTitleTypes);

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

    % sort the rows based on the isCell and ROIindex
    suite2pTable = sortrows(suite2pTable,{'isCell','ROIindexS2P'},{'descend','ascend'});
    save(fullfile(planePath, 'suite2pTable.mat'),'suite2pTable');

    suite2pImage = zeros(512,512);
    for i = 1 : size(suite2pTable, 1)
        if suite2pTable.isCell(i) == 1  % only get the isCell spatial footprint
            suite2pImage = suite2pImage + suite2pTable.cellRegInputS2P{i,1};
        end
    end
    figure; 
    imshow(suite2pImage, []);
    title(['Spatial footprint of suite2p (isCell only):', subDirs(k).name], 'FontSize', 16);


    % get the iscell 3D matrix for cellReg input, save to cellRegInput.mat

    isCellNum = sum(suite2pTable.isCell);
    suite2pInput = zeros(512,512,isCellNum);

    for i = 1 : isCellNum
        suite2pInput(:,:,i) = suite2pTable.cellRegInputS2P{i};
    end

    cellRegInput = permute(suite2pInput,[3 2 1]);

    %use the binary mask to run the cellReg 24.10.31
    cellRegInput = double(cellRegInput ~= 0);

    saveStr = fullfile(planePath, 'cellRegInput.mat');
    save(saveStr,"cellRegInput");

    disp(['Processing completed for folder: ', planePath]);
end
close(wb)