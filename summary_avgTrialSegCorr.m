%% start analysis
clearvars -except ROISegTraceTable_*
%ask for the important parameter settings
dlgTitle = 'Set the threshold parameters';
promptPara = {'Set the CC threshold','Set the spatial distance threshold'};
fieldSize = [1 100;1 100];
defInput = {'0.5','20'};
thresholds = inputdlg(promptPara,dlgTitle,fieldSize,defInput);

ccThreshold = str2double(thresholds{1});
spatialThreshold = str2double(thresholds{2});

contrastLevel = 0.001; % 
behavioralResult = 4; %Hit trial only
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
        tempHitTrace1{j} = interp1(tempHitTrace1{j}(:,1),tempHitTrace1{j}(:,2),commonTimeAxis,'linear','extrap'); %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace1{:});                            %catenate the traces of 100% hit trials
 %calculate the average of each cell in the tempHitTrace1
    ROINum_M3_100HitTrace{i,1} = mean(hitTrace)';
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

%% 
% Analysis for M9
ROINum_M9 = height(ROISegTraceTable_M9);
ROINum_M9_100HitTrace = cell(ROINum_M9,1); %store the hit trace of each ROI in M9

for i = 1:ROINum_M9
    tempTraceTable = ROISegTraceTable_M9.Segmented_trace{i,1};   %get the segmented trace table of current ROI
    tempHitTrace = tempTraceTable.traces(tempTraceTable.trialResult == behavioralResult & tempTraceTable.trialContrast == contrastLevel,:);  %100% hit trials
    tempHitTrace = cellfun(@(x) x(:,2:3),tempHitTrace,'UniformOutput',false);  %get the 3rd column of each cell
    %interpolate the neuronal activity on a common time axis
    for j = 1:length(tempHitTrace)
        tempHitTrace{j} = interp1(tempHitTrace{j}(:,1),tempHitTrace{j}(:,2),commonTimeAxis,'linear','extrap'); %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace{:});                            %catenate the traces of 100% hit trials
    ROINum_M9_100HitTrace{i,1} = mean(hitTrace)';
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

%% Analysis for M17
ROINum_M17 = height(ROISegTraceTable_M17);
ROINum_M17_100HitTrace = cell(ROINum_M17,1); %store the hit trace of each ROI in M17

for i = 1:ROINum_M17
    tempTraceTable = ROISegTraceTable_M17.Segmented_trace{i,1};   %get the segmented trace table of current ROI
    tempHitTrace = tempTraceTable.traces(tempTraceTable.trialResult == behavioralResult & tempTraceTable.trialContrast == contrastLevel,:);  %100% hit trials
    tempHitTrace = cellfun(@(x) x(:,2:3),tempHitTrace,'UniformOutput',false);  %get the 3rd column of each cell
    %interpolate the neuronal activity on a common time axis
    for j = 1:length(tempHitTrace)
        tempHitTrace{j} = interp1(tempHitTrace{j}(:,1),tempHitTrace{j}(:,2),commonTimeAxis,'linear','extrap'); %interpolate the trace
    end
    
    hitTrace = vertcat(tempHitTrace{:});                            %catenate the traces of 100% hit trials
    ROINum_M17_100HitTrace{i,1} = mean(hitTrace)';
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
