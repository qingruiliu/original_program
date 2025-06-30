% check the traces appeared across 3 sessions and checked the preference of individual ROIs

%load the SegTraceTable_session1 .mat file
disp('---select the SegTraceTable_session1.mat file---');
[segTraceFile1,segTracePath1] = uigetfile('*_session1.mat');
cd(segTracePath1);
load(segTraceFile1);
segTraceTable1 = ROISegTraceTable;
segTraceIdx1 = segTraceTable1.ROIIndex; %get the ROI index of session 1

%load the SegTraceTable_session2 .mat file
disp('---select the SegTraceTable_session2.mat file---');
[segTraceFile2,segTracePath2] = uigetfile('*_session2.mat');
cd(segTracePath2);
load(segTraceFile2);
segTraceTable2 = ROISegTraceTable;
segTraceIdx2 = segTraceTable2.ROIIndex; %get the ROI index of session 2

%load the SegTraceTable_session3 .mat file
disp('---select the SegTraceTable_session3.mat file---');
[segTraceFile3,segTracePath3] = uigetfile('*_session3.mat');
cd(segTracePath3);
load(segTraceFile3);
segTraceTable3 = ROISegTraceTable;
segTraceIdx3 = segTraceTable3.ROIIndex; %get the ROI index of session 3

%% find the common ROIs across 3 sessions
commonROIs = intersect(intersect(segTraceIdx1, segTraceIdx2), segTraceIdx3);
disp(['Number of common ROIs across 3 sessions: ', num2str(length(commonROIs))]);

%build a table to store the common ROIs and their traces
commonROITable = table('Size', [length(commonROIs), 4], ...
                        'VariableTypes', {'double', 'cell', 'cell', 'cell'}, ...
                        'VariableNames', {'CommonROIIndex', 'TraceSession1', 'TraceSession2', 'TraceSession3'});

% Fill the table with the common ROIs and their traces
for i = 1:length(commonROIs)
    commonROITable.CommonROIIndex(i) = commonROIs(i);
    commonROITable.TraceSession1{i} = segTraceTable1.Segmented_trace{segTraceIdx1 == commonROIs(i)};
    commonROITable.TraceSession2{i} = segTraceTable2.Segmented_trace{segTraceIdx2 == commonROIs(i)};
    commonROITable.TraceSession3{i} = segTraceTable3.Segmented_trace{segTraceIdx3 == commonROIs(i)};
    
end

%plot the 100% hit trials traces across 3 sessions
for i = 1:length(commonROIs)
    %adjust the height of the figure to fit 3 subplots
    set(0, 'DefaultFigurePosition', [100, 100, 800, 1200]); %set the default figure position
    % Create a new figure for each ROI  
    figure;

    % Define common interpolated time axis
    interpTime = -0.5:0.1:8;

    % Session 1
    subplot(3,1,1);
    hold on;
    currentROITraceTable1 = commonROITable.TraceSession1{i};
    hit100session1Idx = []; % Initialize
    if isfield(currentROITraceTable1, 'trialResult') && isfield(currentROITraceTable1, 'trialContrast')
        hit100session1Idx = find(currentROITraceTable1.trialResult == 1 & currentROITraceTable1.trialContrast == 1);
    else
        disp(['Warning: trialResult/trialContrast fields missing for ROI ', num2str(commonROIs(i)), ' in Session 1.']);
    end
    
    interpolatedTraces1 = [];
    if ~isempty(hit100session1Idx) && isfield(currentROITraceTable1, 'traces') && iscell(currentROITraceTable1.traces) && ~isempty(currentROITraceTable1.traces)
        valid_indices1 = hit100session1Idx(hit100session1Idx <= numel(currentROITraceTable1.traces));
        if ~isempty(valid_indices1)
            hit100session1Trace_subset = currentROITraceTable1.traces(valid_indices1);
            for j = 1:length(hit100session1Trace_subset)
                if ~isempty(hit100session1Trace_subset{j}) && size(hit100session1Trace_subset{j}, 2) >= 3 % Ensure trace has at least time, z-score
                    traceTime = hit100session1Trace_subset{j}(:,2);
                    traceZscore = hit100session1Trace_subset{j}(:,3);
                    [traceTime, uniqueIdx] = unique(traceTime); % Ensure time is monotonically increasing
                    traceZscore = traceZscore(uniqueIdx);
                    if length(traceTime) > 1 % interp1 requires at least 2 data points
                        interpZscore = interp1(traceTime, traceZscore, interpTime, 'linear', NaN);
                        plot(interpTime, interpZscore, 'Color', [0.5 0.5 0.5]);
                        interpolatedTraces1 = [interpolatedTraces1; interpZscore];
                    end
                end
            end
        end
    end
    if ~isempty(interpolatedTraces1)
        meanTrace1 = nanmean(interpolatedTraces1, 1);
        plot(interpTime, meanTrace1, 'r', 'LineWidth', 1.5);
    end
    title(['Session 1 - ROI ', num2str(commonROIs(i))]);
    xlim([-0.5 8]);
    ylim([-0.5 5]);
    xline(1,'--','Color',[0 0 0],'LineWidth',1);
    xline(5,'--','Color',[0 0 0],'LineWidth',1);
    xregion(0,1,'FaceAlpha',0.2);
    set(gca,'TickLength',[0 0]);
    set(gca,'FontSize',14)
    xlabel('Time (s)');
    ylabel('z-score');

    % Session 2
    subplot(3,1,2);
    hold on;
    currentROITraceTable2 = commonROITable.TraceSession2{i};
    hit100session2Idx = []; % Initialize
    if isfield(currentROITraceTable2, 'trialResult') && isfield(currentROITraceTable2, 'trialContrast')
        hit100session2Idx = find(currentROITraceTable2.trialResult == 1 & currentROITraceTable2.trialContrast == 1);
    else
        disp(['Warning: trialResult/trialContrast fields missing for ROI ', num2str(commonROIs(i)), ' in Session 2.']);
    end

    interpolatedTraces2 = [];
    if ~isempty(hit100session2Idx) && isfield(currentROITraceTable2, 'traces') && iscell(currentROITraceTable2.traces) && ~isempty(currentROITraceTable2.traces)
        valid_indices2 = hit100session2Idx(hit100session2Idx <= numel(currentROITraceTable2.traces));
        if ~isempty(valid_indices2)
            hit100session2Trace_subset = currentROITraceTable2.traces(valid_indices2);
            for j = 1:length(hit100session2Trace_subset)
                 if ~isempty(hit100session2Trace_subset{j}) && size(hit100session2Trace_subset{j}, 2) >= 3
                    traceTime = hit100session2Trace_subset{j}(:,2);
                    traceZscore = hit100session2Trace_subset{j}(:,3);
                    [traceTime, uniqueIdx] = unique(traceTime);
                    traceZscore = traceZscore(uniqueIdx);
                    if length(traceTime) > 1
                        interpZscore = interp1(traceTime, traceZscore, interpTime, 'linear', NaN);
                        plot(interpTime, interpZscore, 'Color', [0.5 0.5 0.5]);
                        interpolatedTraces2 = [interpolatedTraces2; interpZscore];
                    end
                end
            end
        end
    end
    if ~isempty(interpolatedTraces2)
        meanTrace2 = nanmean(interpolatedTraces2, 1);
        plot(interpTime, meanTrace2, 'r', 'LineWidth', 1.5);
    end
    title(['Session 2 - ROI ', num2str(commonROIs(i))]);
    xlim([-0.5 8]);
    ylim([-0.5 5]);
    xline(1,'--','Color',[0 0 0],'LineWidth',1);
    xline(5,'--','Color',[0 0 0],'LineWidth',1);
    xregion(0,1,'FaceAlpha',0.2);
    set(gca,'TickLength',[0 0]);
    set(gca,'FontSize',14);   
    xlabel('Time (s)');
    ylabel('z-score');

    % Session 3
    subplot(3,1,3);
    hold on;
    currentROITraceTable3 = commonROITable.TraceSession3{i};
    hit100session3Idx = []; % Initialize
    if isfield(currentROITraceTable3, 'trialResult') && isfield(currentROITraceTable3, 'trialContrast')
        hit100session3Idx = find(currentROITraceTable3.trialResult == 1 & currentROITraceTable3.trialContrast == 1);
    else
        disp(['Warning: trialResult/trialContrast fields missing for ROI ', num2str(commonROIs(i)), ' in Session 3.']);
    end
    
    interpolatedTraces3 = [];
    if ~isempty(hit100session3Idx) && isfield(currentROITraceTable3, 'traces') && iscell(currentROITraceTable3.traces) && ~isempty(currentROITraceTable3.traces)
        valid_indices3 = hit100session3Idx(hit100session3Idx <= numel(currentROITraceTable3.traces));
        if ~isempty(valid_indices3)
            hit100session3Trace_subset = currentROITraceTable3.traces(valid_indices3);
            for j = 1:length(hit100session3Trace_subset)
                if ~isempty(hit100session3Trace_subset{j}) && size(hit100session3Trace_subset{j}, 2) >= 3
                    traceTime = hit100session3Trace_subset{j}(:,2);
                    traceZscore = hit100session3Trace_subset{j}(:,3);
                    [traceTime, uniqueIdx] = unique(traceTime);
                    traceZscore = traceZscore(uniqueIdx);
                    if length(traceTime) > 1
                        interpZscore = interp1(traceTime, traceZscore, interpTime, 'linear', NaN);
                        plot(interpTime, interpZscore, 'Color', [0.5 0.5 0.5]);
                        interpolatedTraces3 = [interpolatedTraces3; interpZscore];
                    end
                end
            end
        end
    end
    if ~isempty(interpolatedTraces3)
        meanTrace3 = nanmean(interpolatedTraces3, 1);
        plot(interpTime, meanTrace3, 'r', 'LineWidth', 1.5);
    end
    title(['Session 3 - ROI ', num2str(commonROIs(i))]);
    xlim([-0.5 8]);
    ylim([-0.5 5]);
    xline(1,'--','Color',[0 0 0],'LineWidth',1);
    xline(5,'--','Color',[0 0 0],'LineWidth',1);
    xregion(0,1,'FaceAlpha',0.2);
    set(gca,'TickLength',[0 0]);
    set(gca,'FontSize',14);
    xlabel('Time (s)');
    ylabel('z-score');

    hold off;

    % dialog box to continue or stop
    answer = questdlg('Do you want to continue to the next ROI?', ...
                      'Continue?', ...
                      'Yes', 'No', 'Yes');
    if strcmp(answer, 'No')
        break; % exit the loop if user chooses 'No'
    end
    % Save the figure
    saveStr = append('ROI_', num2str(commonROIs(i)), '_traces.png');
    saveas(gcf, saveStr);
    close(gcf); % close the figure after saving
end

%%  plot the neuronal activity in individual trials based on different result (10% and 100%)
ROINum = size(neuronMatricesOnResult,1);
hitTrialHigh = length(find(cell2mat(neuronMatricesOnResult(:,2)) == 1 & cell2mat(neuronMatricesOnResult(:,3)) > 0.1));
interpTime = -0.5:0.1:8;

for i = 1 : ROINum
    currentROITrials = neuronMatricesOnResult(:,i+3);
    figure;
    titleStr = append('ROI #',num2str(i),'   100% Hit trials');
    title(titleStr)
    interpTrace = [];
    hold on
    for j = 1 : hitTrialHigh
        traceTime = currentROITrials{j}(:,2);
        traceSpk = currentROITrials{j}(:,3);
        interpZscore = interp1(traceTime,traceSpk,interpTime,'linear','extrap');
        % Corrected plot command for individual trials:
        % Use interpTime and interpZscore directly
        % Changed color to a slightly darker gray and added LineWidth
        plot(interpTime, interpZscore, 'Color', [0.6,0.6,0.6], 'LineWidth', 0.5);
        interpTrace = [interpTrace;interpZscore];
    end
    meanTrace = mean(interpTrace,1);
    plot(interpTime,meanTrace,'r','LineWidth',1.5);
    xlim([-0.5 8]);
    ylim([-0.2 5]);
    xline(1,'--','Color',[0 0 0],'LineWidth',1);
    xline(5,'--','Color',[0 0 0],'LineWidth',1);
    xregion(0,1,'FaceAlpha',0.2);
    set(gca,'TickLength',[0 0]);
    set(gca,'FontSize',14)
    xlabel('Time (s)');
    ylabel('Deconvoluted F_F');
     % dialog box to continue or stop
    answer = questdlg('Do you want to continue to the next ROI?', ...
                      'Continue?', ...
                      'Yes', 'No', 'Yes');
    if strcmp(answer, 'No')
        break; % exit the loop if user chooses 'No'
    end
    % Save the figure
    saveStr = append('ROI_', num2str(commonROIs(i)), '_traces.png');
    saveas(gcf, saveStr);
    close(gcf); % close the figure after saving
end

