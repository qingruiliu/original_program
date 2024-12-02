%% load the ROISegTraceTable from 3 animals
% %M3
disp('----------choose the ROISegTraceTable.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M3 = ROISegTraceTable;
% 
% %M9
disp('----------choose the ROISegTraceTable.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M9 = ROISegTraceTable;

% %M17
disp('----------choose the ROISegTraceTable.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M17 = ROISegTraceTable;

%% start analysis
clearvars -except ROISegTraceTable_*
%ask for the important parameter settings
dlgTitle = 'Set the threshold parameters';
promptPara = {'Set the CC threshold',...
              'Set the spatial distance threshold',...
              'Set the contrast level to process (%)',...
              'Set the result flag to process (1:Hit, 2:Miss, 3:FA, 4:CR)'};
fieldSize = [1 100;1 100;1 100;1 100];
defInput = {'0.5','20','100','1'};
thresholds = inputdlg(promptPara,dlgTitle,fieldSize,defInput);

ccThreshold = str2double(thresholds{1});
spatialThreshold = str2double(thresholds{2});

contrastLevel = str2double(thresholds{3})/100;
behavioralResult = str2double(thresholds{4});
switch behavioralResult
    case 1  % 1: Hit
        plotColor = [0.25 0.8 0.25];
    case 2  % 2: Miss
        plotColor = [1 0.54 0.1];
    case 3  % 3: False Alarm
        plotColor = [0.83 0.14 0.14];
    case 4  % 4: Correct Reject
        plotColor = [0.27 0.25 0.8];
end

%clearvars -except ROISegTraceTable_*
%% next step: contaneate the traces in 100% hit trials 
% (and limit the trial condition in the trial before the hit trial)
commonTimeAxis = -1:0.2:9;                 %common time axis for interpolation
%M3
ROINum_M3 = height(ROISegTraceTable_M3);
ROINum_M3_100HitTrace = cell(ROINum_M3,1); %store the hit trace of each ROI in M3

for i = 1:ROINum_M3
    tempTraceTable = ROISegTraceTable_M3.Segmented_trace{i,1};   %get the segmented trace table of current ROI
    tempHitTrace = tempTraceTable.traces(tempTraceTable.trialResult == behavioralResult & tempTraceTable.trialContrast == contrastLevel,:);  %100% hit trials
    tempHitTrace1 = cellfun(@(x) x(:,2:3),tempHitTrace,'UniformOutput',false);  %get the 3rd column of each cell
    %interpolate the neuronal activity on a common time axis
    for j = 1:length(tempHitTrace1)
        tempHitTrace1{j} = interp1(tempHitTrace1{j}(:,1),tempHitTrace1{j}(:,2),commonTimeAxis,'linear','extrap')'; %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace1{:});                            %catenate the traces of 100% hit trials
    ROINum_M3_100HitTrace{i,1} = hitTrace;
end
% create the index mat for all possible pairs(pairs_M3)
[p,q] = ndgrid(1:ROINum_M3,1:ROINum_M3); %get the index of all possible pairs
pairs_M3 = [p(:),q(:)];                  %get all possible pairs
pairs_M3 = pairs_M3(pairs_M3(:,1) < pairs_M3(:,2),:); %remove the repeated pairs, used for the index of ROI pairs

hitCorrMat_M3 = ROINum_M3_100HitTrace(pairs_M3); %trace pairs
transCenterCell_M3 = num2cell(ROISegTraceTable_M3.TransformedCenter,2); %transformed center of each ROI
hitCenterMat_M3 = transCenterCell_M3(pairs_M3); %spatial distance pairs

hitCorrAndDist_M3 = zeros(size(pairs_M3,1),3); %preset the matrix to store the correlation coefficient and spatial distance of each pair

% calculate the correlation coefficient of all possible pairs
for i = 1:size(pairs_M3,1)
    hitCorrAndDist_M3(i,1) = corr(hitCorrMat_M3{i,1},hitCorrMat_M3{i,2}); %1st column: correlation coefficient
    hitCorrAndDist_M3(i,2) = norm(hitCenterMat_M3{i,1}(1:2) - hitCenterMat_M3{i,2}(1:2)); %2nd column: tangential distance
    hitCorrAndDist_M3(i,3) = norm(hitCenterMat_M3{i,1} - hitCenterMat_M3{i,2}); %3rd column: XYZ spatial distance
end

% figure;
% hold on
% xlim([0 50]);
% ylim([-0.05 0.2]);
% scatter(hitCorrAndDist_M3(:,2),hitCorrAndDist_M3(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
% title('CC and tangential distance in 100% hit trials'); %add the title
% set(gca,'FontSize',16)
% % plot the average correlation coefficient of ROIs with the different distance bin
% distanceBin = 0:5:max(hitCorrAndDist_M3(:,2));       %create the distance bin
% corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
% for i = 1 : length(distanceBin)-1
%     tempIdx = find(hitCorrAndDist_M3(:,2) >= distanceBin(i) & hitCorrAndDist_M3(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
%     corrBin(i) = mean(hitCorrAndDist_M3(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
% end
% %plot the corrbin at the middle of the distance bin
% plot(distanceBin(1:end-1)+2.5,corrBin,'g','LineWidth',4);  %plot the average correlation coefficient of the distance bin
% %plot(distanceBin(2:end)/2,corrBin,'k','LineWidth',2);  %plot the average correlation coefficient of the distance bin
% xlabel('Distance between the paired ROI centroids');     %add the x-axis label
% ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
% set(gca,'FontSize',16)

%
% Analysis for M9
ROINum_M9 = height(ROISegTraceTable_M9);
ROINum_M9_100HitTrace = cell(ROINum_M9,1); %store the hit trace of each ROI in M9

for i = 1:ROINum_M9
    tempTraceTable = ROISegTraceTable_M9.Segmented_trace{i,1};   %get the segmented trace table of current ROI
    tempHitTrace = tempTraceTable.traces(tempTraceTable.trialResult == behavioralResult & tempTraceTable.trialContrast == contrastLevel,:);  %100% hit trials
    tempHitTrace = cellfun(@(x) x(:,2:3),tempHitTrace,'UniformOutput',false);  %get the 3rd column of each cell
    %interpolate the neuronal activity on a common time axis
    for j = 1:length(tempHitTrace)
        tempHitTrace{j} = interp1(tempHitTrace{j}(:,1),tempHitTrace{j}(:,2),commonTimeAxis,'linear','extrap')'; %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace{:});                            %catenate the traces of 100% hit trials
    ROINum_M9_100HitTrace{i,1} = hitTrace;
end
% create the index mat for all possible pairs(pairs_M9)
[p,q] = ndgrid(1:ROINum_M9,1:ROINum_M9); %get the index of all possible pairs
pairs_M9 = [p(:),q(:)];                  %get all possible pairs
pairs_M9 = pairs_M9(pairs_M9(:,1) < pairs_M9(:,2),:); %remove the repeated pairs, used for the index of ROI pairs

hitCorrMat_M9 = ROINum_M9_100HitTrace(pairs_M9); %trace pairs
transCenterCell_M9 = num2cell(ROISegTraceTable_M9.TransformedCenter,2); %transformed center of each ROI
hitCenterMat_M9 = transCenterCell_M9(pairs_M9); %spatial distance pairs

hitCorrAndDist_M9 = zeros(size(pairs_M9,1),3); %preset the matrix to store the correlation coefficient and spatial distance of each pair

% calculate the correlation coefficient of all possible pairs
for i = 1:size(pairs_M9,1)
    hitCorrAndDist_M9(i,1) = corr(hitCorrMat_M9{i,1},hitCorrMat_M9{i,2}); %1st column: correlation coefficient
    hitCorrAndDist_M9(i,2) = norm(hitCenterMat_M9{i,1}(1:2) - hitCenterMat_M9{i,2}(1:2)); %2nd column: tangential distance
    hitCorrAndDist_M9(i,3) = norm(hitCenterMat_M9{i,1} - hitCenterMat_M9{i,2}); %3rd column: XYZ spatial distance
end

% Analysis for M17
ROINum_M17 = height(ROISegTraceTable_M17);
ROINum_M17_100HitTrace = cell(ROINum_M17,1); %store the hit trace of each ROI in M17

for i = 1:ROINum_M17
    tempTraceTable = ROISegTraceTable_M17.Segmented_trace{i,1};   %get the segmented trace table of current ROI
    tempHitTrace = tempTraceTable.traces(tempTraceTable.trialResult == behavioralResult & tempTraceTable.trialContrast == contrastLevel,:);  %100% hit trials
    tempHitTrace = cellfun(@(x) x(:,2:3),tempHitTrace,'UniformOutput',false);  %get the 3rd column of each cell
    %interpolate the neuronal activity on a common time axis
    for j = 1:length(tempHitTrace)
        tempHitTrace{j} = interp1(tempHitTrace{j}(:,1),tempHitTrace{j}(:,2),commonTimeAxis,'linear','extrap')'; %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace{:});                            %catenate the traces of 100% hit trials
    ROINum_M17_100HitTrace{i,1} = hitTrace;
end
% create the index mat for all possible pairs(pairs_M17)
[p,q] = ndgrid(1:ROINum_M17,1:ROINum_M17); %get the index of all possible pairs
pairs_M17 = [p(:),q(:)];                  %get all possible pairs
pairs_M17 = pairs_M17(pairs_M17(:,1) < pairs_M17(:,2),:); %remove the repeated pairs, used for the index of ROI pairs

hitCorrMat_M17 = ROINum_M17_100HitTrace(pairs_M17); %trace pairs
transCenterCell_M17 = num2cell(ROISegTraceTable_M17.TransformedCenter,2); %transformed center of each ROI
hitCenterMat_M17 = transCenterCell_M17(pairs_M17); %spatial distance pairs

hitCorrAndDist_M17 = zeros(size(pairs_M17,1),3); %preset the matrix to store the correlation coefficient and spatial distance of each pair

% calculate the correlation coefficient of all possible pairs
for i = 1:size(pairs_M17,1)
    hitCorrAndDist_M17(i,1) = corr(hitCorrMat_M17{i,1},hitCorrMat_M17{i,2}); %1st column: correlation coefficient
    hitCorrAndDist_M17(i,2) = norm(hitCenterMat_M17{i,1}(1:2) - hitCenterMat_M17{i,2}(1:2)); %2nd column: tangential distance
    hitCorrAndDist_M17(i,3) = norm(hitCenterMat_M17{i,1} - hitCenterMat_M17{i,2}); %3rd column: XYZ spatial distance
end

%% summary the results of 3 animals

hitCorrAndDist_all = [hitCorrAndDist_M3;hitCorrAndDist_M9;hitCorrAndDist_M17]; %combine the results of 3 animals    

% filter the pairs by the threshold
hitCorrAndDist_all = hitCorrAndDist_all(hitCorrAndDist_all(:,3) > spatialThreshold,:); %filter the pairs by spatial distance
hitCorrAndDist_M3 = hitCorrAndDist_all(hitCorrAndDist_all(:,1) <= ccThreshold,:); %filter the pairs by correlation coefficient

figure;
hold on
xlim([0 70]);
xticks(0:20:60);
ylim([-0.03 0.25]);
yticks(0:0.05:0.25);

%scatter(hitCorrAndDist_all(:,2),hitCorrAndDist_all(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('CC and tangential Distance of ROIs in 3 mice'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(hitCorrAndDist_all(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(hitCorrAndDist_all(:,2) >= distanceBin(i) & hitCorrAndDist_all(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(hitCorrAndDist_all(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(1:end-1)+2.5,corrBin,'Color',plotColor,'LineWidth',4);  %plot the average correlation coefficient of the distance bin
xlabel('Tangential distance (µm)');     %add the x-axis label
ylabel('Average correlation');          %add the y-axis label

set(gca,'FontSize',16)

%% plot the surrogate data on the same figure with actual data
surrogateTimes = 1000;
hitSurroMat = zeros(length(pairs_M3) + length(pairs_M9) + length(pairs_M17), 3, surrogateTimes);

wb = waitbar(0, 'Generating the surrogate data...');
parfor i = 1:surrogateTimes
    %M3
    surroPairs_M3_1 = randi(ROINum_M3,size(pairs_M3)); %create the first surrogate pair index
    surroPairs_M3_2 = randi(ROINum_M3,size(pairs_M3)); %create the second surrogate pair index
    
    hitSurroTraceMat_M3 = ROINum_M3_100HitTrace(surroPairs_M3_1); %surrogate trace pairs
    hitSurroCenterMat_M3 = transCenterCell_M3(surroPairs_M3_2); %surrogate spatial distance pairs

    %calculate the correlation coefficient between the every pair between column1 and column2
    hitCorrSurro_M3 = zeros(length(pairs_M3), 3);
    for j = 1:length(pairs_M3)
        hitCorrSurro_M3(j, 1) = corr(hitSurroTraceMat_M3{j,1}, hitSurroTraceMat_M3{j,2}); %calculate the correlation coefficient between the trace of the pair
        hitCorrSurro_M3(j, 2) = norm(hitSurroCenterMat_M3{j,1}(1:2) - hitSurroCenterMat_M3{j,2}(1:2)); %calculate the tangential distance between the pair
        hitCorrSurro_M3(j, 3) = norm(hitSurroCenterMat_M3{j,1} - hitSurroCenterMat_M3{j,2}); %calculate the untransformed XYZ distance between the pair
    end

    %M9
    surroPairs_M9_1 = randi(ROINum_M9,size(pairs_M9)); %create the first surrogate pair index
    surroPairs_M9_2 = randi(ROINum_M9,size(pairs_M9)); %create the second surrogate pair index
    
    hitSurroTraceMat_M9 = ROINum_M9_100HitTrace(surroPairs_M9_1); %surrogate trace pairs
    hitSurroCenterMat_M9 = transCenterCell_M9(surroPairs_M9_2); %surrogate spatial distance pairs

    %calculate the correlation coefficient between the every pair between column1 and column2
    hitCorrSurro_M9 = zeros(length(pairs_M9), 3);
    for j = 1:length(pairs_M9)
        hitCorrSurro_M9(j, 1) = corr(hitSurroTraceMat_M9{j,1}, hitSurroTraceMat_M9{j,2}); %calculate the correlation coefficient between the trace of the pair
        hitCorrSurro_M9(j, 2) = norm(hitSurroCenterMat_M9{j,1}(1:2) - hitSurroCenterMat_M9{j,2}(1:2)); %calculate the tangential distance between the pair
        hitCorrSurro_M9(j, 3) = norm(hitSurroCenterMat_M9{j,1} - hitSurroCenterMat_M9{j,2}); %calculate the untransformed XYZ distance between the pair
    end

    %M17
    surroPairs_M17_1 = randi(ROINum_M17,size(pairs_M17)); %create the first surrogate pair index
    surroPairs_M17_2 = randi(ROINum_M17,size(pairs_M17)); %create the second surrogate pair index
    
    hitSurroTraceMat_M17 = ROINum_M17_100HitTrace(surroPairs_M17_1); %surrogate trace pairs
    hitSurroCenterMat_M17 = transCenterCell_M17(surroPairs_M17_2); %surrogate spatial distance pairs

    %calculate the correlation coefficient between the every pair between column1 and column2
    hitCorrSurro_M17 = zeros(length(pairs_M17), 3);
    for j = 1:length(pairs_M17)
        hitCorrSurro_M17(j, 1) = corr(hitSurroTraceMat_M17{j,1}, hitSurroTraceMat_M17{j,2}); %calculate the correlation coefficient between the trace of the pair
        hitCorrSurro_M17(j, 2) = norm(hitSurroCenterMat_M17{j,1}(1:2) - hitSurroCenterMat_M17{j,2}(1:2)); %calculate the tangential distance between the pair
        hitCorrSurro_M17(j, 3) = norm(hitSurroCenterMat_M17{j,1} - hitSurroCenterMat_M17{j,2}); %calculate the untransformed XYZ distance between the pair
    end

    hitSurroMat(:, :, i) = [hitCorrSurro_M3; hitCorrSurro_M9; hitCorrSurro_M17];
end
close(wb)
%% plot the surrogate data with the actual data figure

surrodistanceBin = 0:5:600;                          %create the distance bin
surrocorrBin = zeros(surrogateTimes,length(surrodistanceBin)-1);       %create a vector to store the average correlation coefficient of each distance bin

wb2 = waitbar(0,'Plotting the surrogate data...');
for i = 1:surrogateTimes
    waitbar(i/surrogateTimes,wb2);

    tempSurroMat = hitSurroMat(:,:,i);
    tempSurroMat = tempSurroMat(tempSurroMat(:,3) > spatialThreshold,:);   %exclude the cell pair with xyz distance less than 30
    tempSurroMat =tempSurroMat(tempSurroMat(:,1) <= ccThreshold,:);  %exclude the cell pair with more than 0.5 CC
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
%plot(distanceBin(1:end-1)+2.5,corrBin,'Color',[0.1 0.7 0.1],'LineWidth',4);  %plot the average correlation coefficient of the distance bin
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
histogram(diffSurro10_20,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_<_1_0 _u_m - CC_1_0_-_2_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 200])
hold on
xline(diffRealCorrBin10_20,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro10_30,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_<_1_0 _u_m - CC_2_0_-_3_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 200])
hold on
xline(diffRealCorrBin10_30,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro10_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_<_1_0 _u_m - CC_3_0_-_4_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.1])
ylim([0 200])
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
histogram(diffSurro20_30,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_1_0_-_2_0 _u_m - CC_2_0_-_3_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 200])
hold on
xline(realCorrBin20_30,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro20_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_1_0_-_2_0 _u_m - CC_3_0_-_4_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 200])
hold on
xline(realCorrBin20_40,'-','Color',[0.1 0.7 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)

figure;
histogram(diffSurro30_40,20,'FaceColor',[0.8 0.8 0.8]);
xlabel('CC_2_0_-_3_0 _u_m - CC_3_0_-_4_0 _u_m');
ylabel('Number of Surrogates');
xlim([-0.02 0.05])
ylim([0 200])
hold on
xline(realCorrBin30_40,'-','Color',[0.7 0.1 0.1],'LineWidth',2);
hold off
set(gca,'FontSize',16)