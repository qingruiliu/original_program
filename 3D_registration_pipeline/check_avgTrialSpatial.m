%% load the ROISegTraceTable from 3 animals 
%M3
disp('----------choose the ROISegTraceTable.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M3 = ROISegTraceTable;

%M9
disp('----------choose the ROISegTraceTable.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M9 = ROISegTraceTable;

%M17
disp('----------choose the ROISegTraceTable.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('ROISegTrace*.mat');
cd(path)
load(fileName)
ROISegTraceTable_M17 = ROISegTraceTable;

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

contrastLevel = 1; % 
behavioralResult = 1; %Hit trial only
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
 %calculate the responsive value of each neuron: -0.5~0 s as the basic activity, 0.5~1 s as the responsive activity
    basicIdx = find(commonTimeAxis >= -0.5 & commonTimeAxis < 0);
    responsiveIdx = find(commonTimeAxis >= 0.5 & commonTimeAxis < 1);
    ROINum_M3_100HitTrace{i,2} = mean(ROINum_M3_100HitTrace{i,1}(responsiveIdx)) - mean(ROINum_M3_100HitTrace{i,1}(basicIdx));
end

%% create the index mat for all possible pairs(pairs_M3)
[p,q] = ndgrid(1:ROINum_M3,1:ROINum_M3); %get the index of all possible pairs
pairs_M3 = [p(:),q(:)];                  %get all possible pairs
pairs_M3 = pairs_M3(pairs_M3(:,1) < pairs_M3(:,2),:); %remove the repeated pairs, used for the index of ROI pairs

ROINum_M3_AvgTrace = ROINum_M3_100HitTrace(:,1);
ROINum_M3_responsiveValue = ROINum_M3_100HitTrace(:,2);
hitCorrMat_M3 = ROINum_M3_AvgTrace(pairs_M3); %trace pairs
transCenterCell_M3 = num2cell(ROISegTraceTable_M3.TransformedCenter,2); %transformed center of each ROI
hitCenterMat_M3 = transCenterCell_M3(pairs_M3); %spatial distance pairs
hitCorrMat_M3(:,3:4) = ROINum_M3_responsiveValue(pairs_M3);
hitCorrMat_M3(:,5:6) = hitCenterMat_M3;


hitCorrAndDist_M3 = zeros(size(pairs_M3,1),3); %preset the matrix to store the correlation coefficient and spatial distance of each pair
stimAndRWIdx = find(commonTimeAxis >= 0 & commonTimeAxis < 5)'; %index of the stimulus and reward period
% calculate the correlation coefficient of all possible pairs
for i = 1:size(pairs_M3,1)
    %only calculate the correlation coefficient during 0~5 seconds
    hitCorrAndDist_M3(i,1) = corr(hitCorrMat_M3{i,1}(stimAndRWIdx),hitCorrMat_M3{i,2}(stimAndRWIdx)); %1st column: correlation coefficient
    hitCorrAndDist_M3(i,2) = norm(hitCenterMat_M3{i,1}(1:2) - hitCenterMat_M3{i,2}(1:2)); %2nd column: tangential distance
    hitCorrAndDist_M3(i,3) = norm(hitCenterMat_M3{i,1} - hitCenterMat_M3{i,2}); %3rd column: XYZ spatial distance
end
hitCorrMat_M3(:,7:9) = num2cell(hitCorrAndDist_M3); %store the correlation coefficient and spatial distance of each pair

%chage the hitCorrMat_M3 as a table
hitCorrMat_M3Table = cell2table(hitCorrMat_M3,'VariableNames',{'trace1','trace2','respV1','respV2','center1','center2','corr','tangentialDist','spatialDist'});

%sort the table by ascending tangentialDist and descending correlation coefficient
hitCorrMat_M3Table1 = sortrows(hitCorrMat_M3Table,{'tangentialDist','corr'},{'ascend','descend'});

%remove the spatially close and highly correlated pairs
%hitCorrMat_M3Table1 = hitCorrMat_M3Table1(hitCorrMat_M3Table1.spatialDist > spatialThreshold,:);
%hitCorrMat_M3Table1 = hitCorrMat_M3Table1(hitCorrMat_M3Table1.corr < ccThreshold,:);

%plot all centroids in a centered spatial volume and highlight the pairs with tangential distance < 20 um and correlation coefficient > 0.1
figure
scatter3(hitCorrMat_M3Table1.center1(:,1),hitCorrMat_M3Table1.center1(:,2),hitCorrMat_M3Table1.center1(:,3),50,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.5)
hold on
axis equal
colorbar
grid off
%% plot the pairs with tangential distance < 20 um and correlation coefficient > 0.1
wb = waitbar(0,'Plotting the pairs with tangential distance < 20 um and correlation coefficient > 0.1...');
% Define the colormap
cmap = parula;
hitCorrMat_M3Table2 = hitCorrMat_M3Table1(hitCorrMat_M3Table1.respV1 >= 0.5 | hitCorrMat_M3Table1.respV2 >= 0.5,:);
responsivePairCount = 0;
for i = 1:height(hitCorrMat_M3Table2)
    waitbar(i/height(hitCorrMat_M3Table2))
    if hitCorrMat_M3Table2.tangentialDist(i) < 20 && hitCorrMat_M3Table2.corr(i) > 0.1
        % Get the color from the colormap based on the original correlation value
        color = cmap(round(hitCorrMat_M3Table2.corr(i) * (size(cmap, 1) - 1)) + 1, :);
        
        % Plot the paired centers with the scaled color
        scatter3(hitCorrMat_M3Table2.center1(i,1), hitCorrMat_M3Table2.center1(i,2), hitCorrMat_M3Table2.center1(i,3), 60, 'filled', 'MarkerFaceColor', color)
        scatter3(hitCorrMat_M3Table2.center2(i,1), hitCorrMat_M3Table2.center2(i,2), hitCorrMat_M3Table2.center2(i,3), 60, 'filled', 'MarkerFaceColor', color)
        %draw a line between the paired centers
        plot3([hitCorrMat_M3Table2.center1(i,1),hitCorrMat_M3Table2.center2(i,1)],[hitCorrMat_M3Table2.center1(i,2),hitCorrMat_M3Table2.center2(i,2)],[hitCorrMat_M3Table2.center1(i,3),hitCorrMat_M3Table2.center2(i,3)],'Color',color)
        responsivePairCount = responsivePairCount + 1;
    end
end
close(wb)
title('Avg 100% Hit trials, tangential<20, 0.1<CC<0.5','FontSize',16)
disp('Responsive cell pairs number:')
disp(responsivePairCount)
