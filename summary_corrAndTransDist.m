%% ask for the important parameters setting

dlgTitle = 'Set the threshold parameters';
promptPara = {'\fontsize{20}Set the CC threshold','\fontsize{20}set the spatial distance threshold'};
fieldSize = [1 100;1 100];
defInput = {'0.5','20'};
thresholds = inputdlg(promptPara,dlgTitle,fieldSize,defInput);

ccThreshold = str2double(thresholds{1});
spatialThreshold = str2double(thresholds{2});

%% summarize and analyze the data from 3 mice
disp('----------choose the corrAndDist.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M3 = corrAndDist;
%delete all the odd rows
CorrAndDist_M3 = CorrAndDist_M3(CorrAndDist_M3(:,2)~=0,:);
CorrAndDist_M3(1:2:end,:) = [];
clearvars CorrAndDist
%
disp('----------choose the corrAndDist.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M9 = corrAndDist;
%delete all the odd rows
CorrAndDist_M9 = CorrAndDist_M9(CorrAndDist_M9(:,2)~=0,:);
CorrAndDist_M9(1:2:end,:) = [];
clearvars CorrAndDist
%
disp('----------choose the corrAndDist.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M17 = corrAndDist;
%delete all the odd rows
CorrAndDist_M17 = CorrAndDist_M17(CorrAndDist_M17(:,2)~=0,:);
CorrAndDist_M17(1:2:end,:) = [];
clearvars CorrAndDist

% cat the data
corrAndDist_all = [CorrAndDist_M3; CorrAndDist_M9; CorrAndDist_M17];

%exclude the same ROI pair
corrAndDist_all = corrAndDist_all(corrAndDist_all(:,2) ~= 0,:);

%exclude the cell pair with more than 0.5 CC
corrAndDist_all = corrAndDist_all(corrAndDist_all(:,1) <= ccThreshold,:);

%sort the data based on the second column(untransformed XY distance)
corrAndDist_all = sortrows(corrAndDist_all,2);

%exclude the cell pair with xyz distance less than 30
corrAndDist_all1 = corrAndDist_all(corrAndDist_all(:,3) > spatialThreshold,:);

%% plot the untransformed distance of the cell pairs
figure;
hold on
xlim([0 50]);
ylim([-0.1 0.2]);
%scatter(corrAndDist_all1(:,2),corrAndDist_all1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('CC and untransformed Distance of ROIs in 3 mice'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndDist_all1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndDist_all1(:,2) >= distanceBin(i) & corrAndDist_all1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndDist_all1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
%plot the corrbin at the middle of the distance bin
plot(distanceBin(1:end-1)+2.5,corrBin,'k','LineWidth',4);  %plot the average correlation coefficient of the distance bin
%plot(distanceBin(2:end)/2,corrBin,'k','LineWidth',2);  %plot the average correlation coefficient of the distance bin
xlabel('Distance between the paired ROI centroids');     %add the x-axis label
ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
set(gca,'FontSize',16)

%% summarize and analyze the data from 3 mice

disp('----------choose the corrAndTransDist.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('corrAndTrans*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M3 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M3 = CorrAndTransDist_M3(CorrAndTransDist_M3(:,2)~=0,:);
CorrAndTransDist_M3(1:2:end,:) = [];
clearvars CorrAndTransDist
%
disp('----------choose the corrAndTransDist.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('corrAndTransDist*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M9 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M9 = CorrAndTransDist_M9(CorrAndTransDist_M9(:,2)~=0,:);
CorrAndTransDist_M9(1:2:end,:) = [];
clearvars CorrAndTransDist
%
disp('----------choose the corrAndTransDist.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('corrAndTransDist*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M17 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M17 = CorrAndTransDist_M17(CorrAndTransDist_M17(:,2)~=0,:);
CorrAndTransDist_M17(1:2:end,:) = [];
clearvars CorrAndTransDist

%% cat the data
corrAndTransDist_all = [CorrAndTransDist_M3; CorrAndTransDist_M9; CorrAndTransDist_M17];

%exclude the same ROI pair
corrAndTransDist_all = corrAndTransDist_all(corrAndTransDist_all(:,2) ~= 0,:);

%exclude the cell pair with 1 correlation coefficient
corrAndTransDist_all = corrAndTransDist_all(corrAndTransDist_all(:,1) <= ccThreshold,:);
%sort the data based on the second column(transDist)
corrAndTransDist_all = sortrows(corrAndTransDist_all,2);

%exclude the cell pair with xyz distance less than 30
corrAndTransDist_all1 = corrAndTransDist_all(corrAndTransDist_all(:,3) > spatialThreshold,:);

figure;
hold on
xlim([0 70]);
xticks(0:20:60);
ylim([-0.03 0.12]);
yticks(0:0.05:0.1);

%scatter(corrAndTransDist_all1(:,2),corrAndTransDist_all1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('CC and tangential Distance of ROIs in 3 mice'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndTransDist_all1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndTransDist_all1(:,2) >= distanceBin(i) & corrAndTransDist_all1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndTransDist_all1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(1:end-1)+2.5,corrBin,'Color',[0.1 0.7 0.1],'LineWidth',4);  %plot the average correlation coefficient of the distance bin
xlabel('Tangential distance (µm)');     %add the x-axis label
ylabel('Average correlation');          %add the y-axis label

set(gca,'FontSize',16)

%% load the table of 3 mice and generate the surrogate data 
disp('----------choose the ROICenterAndTraceTable.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('ROICenter*.mat');
cd(path)
load(fileName)
ROICorrIdxTable_M3 = ROICorrIdxTable;
ROINum_M3 = size(ROICorrIdxTable,1);
%
disp('----------choose the ROICenterAndTraceTable.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('ROICenter*.mat');
cd(path)
load(fileName)
ROICorrIdxTable_M9 = ROICorrIdxTable;
ROINum_M9 = size(ROICorrIdxTable,1);
%
disp('----------choose the ROICenterAndTraceTable.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('ROICenter*.mat');
cd(path)
load(fileName)
ROICorrIdxTable_M17 = ROICorrIdxTable;
ROINum_M17 = size(ROICorrIdxTable,1);

%% pre-set the variable for surrogate data
surrogateTimes = 1000;
surroMat = zeros(ROINum_M3*(ROINum_M3-1)/2 + ROINum_M9*(ROINum_M9-1)/2 + ROINum_M17*(ROINum_M17-1)/2, 3, surrogateTimes);

wb = waitbar(0, 'Generating the surrogate data...');
parfor i = 1:surrogateTimes
    %M3
    randIdx_M3 = randi(ROINum_M3, ROINum_M3*(ROINum_M3-1)/2, 4);   %generate the 4-column random index for the 1st animal
    traceSurro_M3 = ROICorrIdxTable_M3.Trace(randIdx_M3(:, 1:2));    %get the trace of the random index
    transCenterSurro_M3_1 = ROICorrIdxTable_M3.TransformedCenter(randIdx_M3(:,3), :);
    transCenterSurro_M3_2 = ROICorrIdxTable_M3.TransformedCenter(randIdx_M3(:,4), :);

    %calculate the correlation coefficient between the every pair between column1 and column2
    corrSurro_M3 = zeros(ROINum_M3*(ROINum_M3-1)/2, 3);
    for j = 1:ROINum_M3*(ROINum_M3-1)/2
        corrSurro_M3(j, 1) = corr(traceSurro_M3{j, 1}, traceSurro_M3{j, 2}); %calculate the correlation coefficient between the trace of the pair
        corrSurro_M3(j, 2) = norm(transCenterSurro_M3_1(j,1:2) - transCenterSurro_M3_2(j,1:2)); %calculate the tangential distance between the pair
        corrSurro_M3(j, 3) = norm(transCenterSurro_M3_1(j,:) - transCenterSurro_M3_2(j,:)); %calculate the untransformed XYZ distance between the pair
    end

    %M9
    randIdx_M9 = randi(ROINum_M9, ROINum_M9*(ROINum_M9-1)/2, 4);   %generate the 4-column random index for the 2nd animal
    traceSurro_M9 = ROICorrIdxTable_M9.Trace(randIdx_M9(:, 1:2));    %get the trace of the random index
    transCenterSurro_M9_1 = ROICorrIdxTable_M9.TransformedCenter(randIdx_M9(:,3), :);
    transCenterSurro_M9_2 = ROICorrIdxTable_M9.TransformedCenter(randIdx_M9(:,4), :);

    %calculate the correlation coefficient between the every pair between column1 and column2
    corrSurro_M9 = zeros(ROINum_M9*(ROINum_M9-1)/2, 3);
    for j = 1:ROINum_M9*(ROINum_M9-1)/2
        corrSurro_M9(j, 1) = corr(traceSurro_M9{j, 1}, traceSurro_M9{j, 2}); %calculate the correlation coefficient between the trace of the pair
        corrSurro_M9(j, 2) = norm(transCenterSurro_M9_1(j,1:2) - transCenterSurro_M9_2(j,1:2)); %calculate the tangential distance between the pair
        corrSurro_M9(j, 3) = norm(transCenterSurro_M9_1(j,:) - transCenterSurro_M9_2(j,:)); %calculate the untransformed XYZ distance between the pair
    end

    %M17
    randIdx_M17 = randi(ROINum_M17, ROINum_M17*(ROINum_M17-1)/2, 4);   %generate the 4-column random index for the 3rd animal
    traceSurro_M17 = ROICorrIdxTable_M17.Trace(randIdx_M17(:, 1:2));    %get the trace of the random index
    transCenterSurro_M17_1 = ROICorrIdxTable_M17.TransformedCenter(randIdx_M17(:,3), :);
    transCenterSurro_M17_2 = ROICorrIdxTable_M17.TransformedCenter(randIdx_M17(:,4), :);

    %calculate the correlation coefficient between the every pair between column1 and column2
    corrSurro_M17 = zeros(ROINum_M17*(ROINum_M17-1)/2, 3);
    for j = 1:ROINum_M17*(ROINum_M17-1)/2
        corrSurro_M17(j, 1) = corr(traceSurro_M17{j, 1}, traceSurro_M17{j, 2}); %calculate the correlation coefficient between the trace of the pair
        corrSurro_M17(j, 2) = norm(transCenterSurro_M17_1(j,1:2) - transCenterSurro_M17_2(j,1:2)); %calculate the tangential distance between the pair
        corrSurro_M17(j, 3) = norm(transCenterSurro_M17_1(j,:) - transCenterSurro_M17_2(j,:)); %calculate the untransformed XYZ distance between the pair
    end

    surroMat(:, :, i) = [corrSurro_M3; corrSurro_M9; corrSurro_M17];
end
close(wb)
% plot the surrogate data with the actual data figure

surrodistanceBin = 0:5:600;                          %create the distance bin
surrocorrBin = zeros(surrogateTimes,length(surrodistanceBin)-1);       %create a vector to store the average correlation coefficient of each distance bin

wb2 = waitbar(0,'Plotting the surrogate data...');
for i = 1:surrogateTimes
    waitbar(i/surrogateTimes,wb2);

    tempSurroMat = surroMat(:,:,i);
    tempSurroMat = tempSurroMat(tempSurroMat(:,3) > spatialThreshold,:);   %exclude the cell pair with xyz distance less than 30
    %tempSurroMat =tempSurroMat(tempSurroMat(:,1) <= ccThreshold,:);  %exclude the cell pair with more than 0.5 CC
    %tempSurroMat(1:2:end,:) = [];                            %delete all the odd rows

    for j = 1 : length(surrodistanceBin)-1
        tempIdx = find(tempSurroMat(:,2) >= surrodistanceBin(j) & tempSurroMat(:,2) < surrodistanceBin(j+1)); %find the index of the correlation coefficient in the distance bin
        surrocorrBin(i,j) = mean(tempSurroMat(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
    end
    plot(surrodistanceBin(1:end-1)+2.5,surrocorrBin(i,:),'Color',[0.5 0.5 0.5 0.2],'LineWidth',2);  %plot the average correlation coefficient of the distance bin
end
close(wb2)

% plot the highest 2.5% and lowest 2.5% of the surrogate data
surrocorrBin1 = sort(surrocorrBin,1,'ascend'); %sort the surrogate data along the column
surroUpperBound = surrocorrBin1(ceil(surrogateTimes*0.975),:); %get the highest 2.5% of the surrogate data
surroLowerBound = surrocorrBin1(ceil(surrogateTimes*0.025),:); %get the lowest 2.5% of the surrogate data

plot(surrodistanceBin(1:end-1)+2.5,surroUpperBound,'k--','LineWidth',1); %plot the highest 2.5% of the surrogate data
plot(surrodistanceBin(1:end-1)+2.5,surroLowerBound,'k--','LineWidth',1); %plot the lowest 2.5% of the surrogate data

%plot the mean of the surrogate data
surroMean = mean(surrocorrBin,1); %calculate the mean of the surrogate data
plot(surrodistanceBin(1:end-1),surroMean,'k','LineWidth',1.5); %plot the mean of the surrogate data
plot(distanceBin(1:end-1)+2.5,corrBin,'Color',[0.1 0.7 0.1],'LineWidth',4);  %plot the average correlation coefficient of the distance bin
hold off

%% make the bin plot of the surrogate distribution and the real data

surroBin10 = surrocorrBin(:,find(distanceBin <= 10)); %get the surrogate data <10 um bins
surroBin20 = surrocorrBin(:,find(distanceBin > 10 & distanceBin <= 20 )); 
surroBin30 = surrocorrBin(:,find(distanceBin > 20 & distanceBin <= 30 )); 
surroBin40 = surrocorrBin(:,find(distanceBin > 30 & distanceBin <= 40 )); 

surroBin10mean = mean(surroBin10,2); %get the mean of the surrogate data <10 um bins
surroBin20mean = mean(surroBin20,2); 
surroBin30mean = mean(surroBin30,2); 
surroBin40mean = mean(surroBin40,2); 

diffSurro10_20 = surroBin10mean - surroBin20mean; %get the difference of the mean of the surrogate data <10 um bins and 10-20 um bins
diffSurro10_30 = surroBin10mean - surroBin30mean;
diffSurro10_40 = surroBin10mean - surroBin40mean;

realCorrBin10 = corrBin(find(distanceBin <= 10)); %get the real data <10 um bins
realCorrBin20 = corrBin(find(distanceBin > 10 & distanceBin <= 20 ));
realCorrBin30 = corrBin(find(distanceBin > 20 & distanceBin <= 30 ));
realCorrBin40 = corrBin(find(distanceBin > 30 & distanceBin <= 40 ));

diffRealCorrBin10_20 = mean(realCorrBin10) - mean(realCorrBin20); %get the difference of the mean of the real data <10 um bins and 10-20 um bins
diffRealCorrBin10_30 = mean(realCorrBin10) - mean(realCorrBin30);
diffRealCorrBin10_40 = mean(realCorrBin10) - mean(realCorrBin40);

%plot the histogram of the difference of the mean of the surrogate data

figure;
histogram(diffSurro10_20,20,'FaceColor',[0.8 0.8 0.8],'LineWidth',1.5);
xlabel('CC_<_1_0 _µ_m - CC_1_0_-_2_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 150])
hold on
xline(diffRealCorrBin10_20,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro10_30,20,'FaceColor',[0.8 0.8 0.8],'LineWidth',1.5);
xlabel('CC_<_1_0 _µ_m - CC_2_0_-_3_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 150])
hold on
xline(diffRealCorrBin10_30,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro10_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_<_1_0 _µ_m - CC_3_0_-_4_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 150])
hold on 
xline(diffRealCorrBin10_40,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

% plot the histogram of 20-30,20-40,30-40
diffSurro20_30 = surroBin20mean - surroBin30mean; %get the difference of the mean of the surrogate data 20-30 um bins and 30-40 um bins
diffSurro20_40 = surroBin20mean - surroBin40mean;
diffSurro30_40 = surroBin30mean - surroBin40mean;

realCorrBin20_30 = mean(realCorrBin20) - mean(realCorrBin30); %get the difference of the mean of the real data 20-30 um bins and 30-40 um bins
realCorrBin20_40 = mean(realCorrBin20) - mean(realCorrBin40);
realCorrBin30_40 = mean(realCorrBin30) - mean(realCorrBin40);

figure;
histogram(diffSurro20_30,20,'FaceColor',[0.8 0.8 0.8],'LineWidth',1.5);
xlabel('CC_1_0_-_2_0 _µ_m - CC_2_0_-_3_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 150])
hold on
xline(realCorrBin20_30,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro20_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_1_0_-_2_0 _µ_m - CC_3_0_-_4_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 150])
hold on
xline(realCorrBin20_40,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro30_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_2_0_-_3_0 _µ_m - CC_3_0_-_4_0 _µ_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 150])
hold on
xline(realCorrBin30_40,'-','Color',[0.7 0.1 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)