%% Comprehensive Behavioral Training Analysis Script
% This script analyzes discrimination training data from stage3_Demo240119.m
% Calculates: Correct Rate, Hit Rate, FA Rate, and Hit Latency across training days
% Can also load pre-processed behavioral_analysis_results.mat for visualization
% Author: Modified for Liu's behavioral analysis
% Date: 2026.02.02

%% 1. Choose data source
choice = questdlg('Select data source:', 'Data Source', ...
    'Raw training files', 'Pre-processed results', 'Raw training files');

if isempty(choice)
    disp('Analysis canceled');
    return;
end

metrics = struct();
path = '';

if strcmp(choice, 'Pre-processed results')
    %% Load pre-processed results
    [fileName, path] = uigetfile('*.mat', 'Select behavioral_analysis_results.mat file');
    
    if isequal(fileName, 0)
        disp('User canceled file selection');
        return;
    end
    
    fullPath = fullfile(path, fileName);
    data = load(fullPath);
    
    if ~isfield(data, 'metrics')
        error('Selected file does not contain "metrics" field. Please select a valid behavioral_analysis_results.mat file.');
    end
    
    metrics = data.metrics;
    disp(['Loaded pre-processed results from: ' fileName]);
    nDays = length(metrics.correctRate);
    disp(['Data contains ' num2str(nDays) ' training day(s)']);
    
else
    %% Load and process raw training session files
    [fileNames, path] = uigetfile('*.mat', 'Select training data files (multiple days)', 'MultiSelect', 'on');
    
    if isequal(fileNames, 0)
        disp('User canceled file selection');
        return;
    end
    
    % Ensure fileNames is a cell array
    if ischar(fileNames)
        fileNames = {fileNames};
    end
    
    nDays = numel(fileNames);
    disp(['Loading ' num2str(nDays) ' training session(s)...']);

    nDays = numel(fileNames);
    disp(['Loading ' num2str(nDays) ' training session(s)...']);
    
    %% 2. Initialize storage variables
    % Structure to store all metrics
    metrics.correctRate = zeros(nDays, 1);
    metrics.hitRate = zeros(nDays, 1);
    metrics.FARate = zeros(nDays, 1);
    metrics.CRRate = zeros(nDays, 1);
    metrics.missRate = zeros(nDays, 1);
    
    % Trial counts
    metrics.nHit = zeros(nDays, 1);
    metrics.nMiss = zeros(nDays, 1);
    metrics.nFA = zeros(nDays, 1);
    metrics.nCR = zeros(nDays, 1);
    metrics.nTotal = zeros(nDays, 1);
    metrics.nTarget = zeros(nDays, 1);
    metrics.nNonTarget = zeros(nDays, 1);
    
    % Hit latency data (stored as cell array for variable trial numbers)
    metrics.hitLatency = cell(nDays, 1);
    metrics.hitLatencyMedian = zeros(nDays, 1);
    metrics.hitLatencyStd = zeros(nDays, 1);
    
    % File names for reference
    metrics.fileNames = fileNames;
    
    %% 3. Process each training session
    wb = waitbar(0, 'Processing training sessions...');
    
    for dayIdx = 1:nDays
        waitbar(dayIdx/nDays, wb, sprintf('Processing day %d/%d...', dayIdx, nDays));
        
        % Load data
        fullPath = fullfile(path, fileNames{dayIdx});
        data = load(fullPath);
        h = data.h;
        
        % Extract result flags from h.data1
        % Column 2: 1=Hit, 2=Miss, 3=FA, 4=CR
        resultFlags = h.data1(:, 2);
        resultFlags(resultFlags == 0) = []; % Remove empty trials
        
        % Count each trial type
        metrics.nHit(dayIdx) = sum(resultFlags == 1);
        metrics.nMiss(dayIdx) = sum(resultFlags == 2);
        metrics.nFA(dayIdx) = sum(resultFlags == 3);
        metrics.nCR(dayIdx) = sum(resultFlags == 4);
        metrics.nTotal(dayIdx) = length(resultFlags);
        
        % Calculate target and non-target trials
        metrics.nTarget(dayIdx) = metrics.nHit(dayIdx) + metrics.nMiss(dayIdx);
        metrics.nNonTarget(dayIdx) = metrics.nFA(dayIdx) + metrics.nCR(dayIdx);
        
        % Calculate rates
        % Correct Rate = (Hit + CR) / Total
        metrics.correctRate(dayIdx) = (metrics.nHit(dayIdx) + metrics.nCR(dayIdx)) / metrics.nTotal(dayIdx);
        
        % Hit Rate = Hit / (Hit + Miss)
        if metrics.nTarget(dayIdx) > 0
            metrics.hitRate(dayIdx) = metrics.nHit(dayIdx) / metrics.nTarget(dayIdx);
        else
            metrics.hitRate(dayIdx) = NaN;
        end
        
        % FA Rate = FA / (FA + CR)
        if metrics.nNonTarget(dayIdx) > 0
            metrics.FARate(dayIdx) = metrics.nFA(dayIdx) / metrics.nNonTarget(dayIdx);
        else
            metrics.FARate(dayIdx) = NaN;
        end
        
        % CR Rate = CR / (FA + CR)
        if metrics.nNonTarget(dayIdx) > 0
            metrics.CRRate(dayIdx) = metrics.nCR(dayIdx) / metrics.nNonTarget(dayIdx);
        else
            metrics.CRRate(dayIdx) = NaN;
        end
        
        % Miss Rate = Miss / (Hit + Miss)
        if metrics.nTarget(dayIdx) > 0
            metrics.missRate(dayIdx) = metrics.nMiss(dayIdx) / metrics.nTarget(dayIdx);
        else
            metrics.missRate(dayIdx) = NaN;
        end
        
        % Extract hit trial latencies
        hitLatencies = [];
        for trialIdx = 1:length(resultFlags)
            if resultFlags(trialIdx) == 1 && trialIdx <= length(h.lickdata)
                % Get first lick latency in hit trial
                if ~isempty(h.lickdata{trialIdx}) && size(h.lickdata{trialIdx}, 2) >= 2
                    latency = h.lickdata{trialIdx}(1, 2);
                    % Filter out abnormal values
                    if latency > 0 && latency < 10
                        hitLatencies = [hitLatencies; latency];
                    end
                end
            end
        end
        
        metrics.hitLatency{dayIdx} = hitLatencies;
        
        if ~isempty(hitLatencies)
            metrics.hitLatencyMedian(dayIdx) = median(hitLatencies);
            metrics.hitLatencyStd(dayIdx) = std(hitLatencies);
        else
            metrics.hitLatencyMedian(dayIdx) = NaN;
            metrics.hitLatencyStd(dayIdx) = NaN;
        end
    end
    
    close(wb);
    disp('Data processing complete!');
end

%% 4. Display summary statistics
fprintf('\n===== Training Summary =====\n');
fprintf('Total training days: %d\n', nDays);

if isfield(metrics, 'nTotal')
    fprintf('Total trials: %d\n', sum(metrics.nTotal));
    fprintf('\n');
    
    for dayIdx = 1:nDays
        if isfield(metrics, 'fileNames')
            fprintf('Day %d (%s):\n', dayIdx, metrics.fileNames{dayIdx});
        else
            fprintf('Day %d:\n', dayIdx);
        end
        fprintf('  Trials: %d (Target: %d, Non-target: %d)\n', ...
            metrics.nTotal(dayIdx), metrics.nTarget(dayIdx), metrics.nNonTarget(dayIdx));
        fprintf('  Hit: %d, Miss: %d, FA: %d, CR: %d\n', ...
            metrics.nHit(dayIdx), metrics.nMiss(dayIdx), metrics.nFA(dayIdx), metrics.nCR(dayIdx));
        fprintf('  Correct Rate: %.2f%%\n', metrics.correctRate(dayIdx) * 100);
        fprintf('  Hit Rate: %.2f%%\n', metrics.hitRate(dayIdx) * 100);
        fprintf('  FA Rate: %.2f%%\n', metrics.FARate(dayIdx) * 100);
        fprintf('  Hit Latency: %.3f ± %.3f s\n', ...
            metrics.hitLatencyMedian(dayIdx), metrics.hitLatencyStd(dayIdx));
        fprintf('\n');
    end
else
    fprintf('\n');
    for dayIdx = 1:nDays
        fprintf('Day %d:\n', dayIdx);
        fprintf('  Correct Rate: %.2f%%\n', metrics.correctRate(dayIdx) * 100);
        fprintf('  Hit Rate: %.2f%%\n', metrics.hitRate(dayIdx) * 100);
        fprintf('  FA Rate: %.2f%%\n', metrics.FARate(dayIdx) * 100);
        fprintf('  Hit Latency: %.3f ± %.3f s\n', ...
            metrics.hitLatencyMedian(dayIdx), metrics.hitLatencyStd(dayIdx));
        fprintf('\n');
    end
end

%% 5. Visualization - Create separate figures for each metric
fprintf('Creating individual figures for each metric...\n');

% Store all figure handles
figHandles = struct();

% Color scheme
colorCorrect = [0.2, 0.7, 0.3];
colorHit = [0.3, 0.5, 0.9];
colorFA = [0.9, 0.3, 0.3];
colorLatency = [0.5, 0.2, 0.7];

days = 1:nDays;

% Figure 1: Correct Rate
figHandles.correctRate = figure('Name', 'Correct Rate', 'Position', [100, 100, 700, 550]);
plot(days, metrics.correctRate * 100, 'o-', 'Color', colorCorrect, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorCorrect);
yline(80, '--k', 'LineWidth', 2, 'Alpha', 0.5);
xlabel('Training Day', 'FontSize', 14);
ylabel('Correct Rate (%)', 'FontSize', 14);
title('Correct Rate Across Days', 'FontSize', 16, 'FontWeight', 'bold');
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on;

% Figure 2: Hit Rate and FA Rate
figHandles.hitFARate = figure('Name', 'Hit Rate and FA Rate', 'Position', [150, 150, 700, 550]);
hold on;
plot(days, metrics.hitRate * 100, 'o-', 'Color', colorHit, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorHit, 'DisplayName', 'Hit Rate');
plot(days, metrics.FARate * 100, 's-', 'Color', colorFA, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorFA, 'DisplayName', 'FA Rate');
xlabel('Training Day', 'FontSize', 14);
ylabel('Rate (%)', 'FontSize', 14);
title('Hit Rate and FA Rate', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
ylim([0, 105]);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
box on;
hold off;

% Figure 3: Trial Counts Stacked Bar (only if detailed data available)
if isfield(metrics, 'nHit')
    figHandles.trialCounts = figure('Name', 'Trial Type Distribution', 'Position', [200, 200, 700, 550]);
    trialCounts = [metrics.nHit, metrics.nCR, metrics.nFA, metrics.nMiss];
    bar(days, trialCounts, 'stacked');
    xlabel('Training Day', 'FontSize', 14);
    ylabel('Number of Trials', 'FontSize', 14);
    title('Trial Type Distribution', 'FontSize', 16, 'FontWeight', 'bold');
    legend({'Hit', 'CR', 'FA', 'Miss'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;
end

% Figure 4: Hit Latency Median with Error Bars
figHandles.hitLatencyBar = figure('Name', 'Hit Latency with Error Bars', 'Position', [250, 250, 700, 550]);
validDays = ~isnan(metrics.hitLatencyMedian);
if any(validDays)
    errorbar(days(validDays), metrics.hitLatencyMedian(validDays), metrics.hitLatencyStd(validDays), 'o-', ...
        'Color', colorLatency, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorLatency);
    xlabel('Training Day', 'FontSize', 14);
    ylabel('Hit Latency (s)', 'FontSize', 14);
    title('Hit Latency (Median ± SD)', 'FontSize', 16, 'FontWeight', 'bold');
    maxLatency = max(metrics.hitLatencyMedian(validDays) + metrics.hitLatencyStd(validDays));
    if ~isempty(maxLatency) && maxLatency > 0
        ylim([0, maxLatency * 1.2]);
    end
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;
end

% Figure 5: Hit Latency Distribution (Scatter with jitter)
if isfield(metrics, 'hitLatency')
    figHandles.hitLatencyDist = figure('Name', 'Hit Latency Distribution', 'Position', [300, 300, 700, 550]);
    hold on;
    grayColor = [0.7, 0.7, 0.7];
    for dayIdx = 1:nDays
        if ~isempty(metrics.hitLatency{dayIdx})
            latencies = metrics.hitLatency{dayIdx};
            x = dayIdx * ones(size(latencies));
            jitter = 0.15 * randn(size(x));
            scatter(x + jitter, latencies, 36, grayColor, 'filled', 'MarkerFaceAlpha', 0.4);
            
            % Plot median line
            medVal = metrics.hitLatencyMedian(dayIdx);
            if ~isnan(medVal)
                line([dayIdx - 0.3, dayIdx + 0.3], [medVal, medVal], 'Color', 'k', 'LineWidth', 3);
            end
        end
    end
    xlabel('Training Day', 'FontSize', 14);
    ylabel('Hit Latency (s)', 'FontSize', 14);
    title('Hit Latency Distribution', 'FontSize', 16, 'FontWeight', 'bold');
    ylim([0, 5]);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    box on;
    hold off;
end

% Figure 6: Performance Matrix (Confusion-like)
if isfield(metrics, 'hitRate') && isfield(metrics, 'missRate')
    figHandles.perfMatrix = figure('Name', 'Performance Matrix', 'Position', [350, 350, 650, 550]);
    % Calculate average rates
    avgHitRate = mean(metrics.hitRate, 'omitnan');
    avgMissRate = mean(metrics.missRate, 'omitnan');
    avgFARate = mean(metrics.FARate, 'omitnan');
    avgCRRate = mean(metrics.CRRate, 'omitnan');
    
    perfMatrix = [avgHitRate, avgMissRate; avgFARate, avgCRRate] * 100;
    imagesc(perfMatrix);
    colormap(flipud(gray));
    colorbar;
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Lick', 'No Lick'}, ...
        'YTick', 1:2, 'YTickLabel', {'Target', 'Non-Target'}, 'FontSize', 12);
    title('Average Performance Matrix (%)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Response', 'FontSize', 14);
    ylabel('Stimulus', 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    
    % Add text annotations
    for i = 1:2
        for j = 1:2
            text(j, i, sprintf('%.1f%%', perfMatrix(i,j)), ...
                'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 16, 'FontWeight', 'bold');
        end
    end
end

fprintf('Created %d individual figures\n', length(fieldnames(figHandles)));

%% 6. Save results
% Prompt user to save the analysis results
answer = questdlg('Do you want to save the analysis results?', 'Save Results', 'Yes', 'No', 'Yes');

if strcmp(answer, 'Yes')
    [saveFileName, savePath] = uiputfile('*.mat', 'Save analysis results as', [path filesep 'behavioral_analysis_results.mat']);
    if ~isequal(saveFileName, 0)
        % Save metrics data
        save(fullfile(savePath, saveFileName), 'metrics');
        disp(['Results saved to: ' fullfile(savePath, saveFileName)]);
        
        % Get base name for files
        [~, baseName, ~] = fileparts(saveFileName);
        
        % Save all individual figures
        figNames = fieldnames(figHandles);
        for fIdx = 1:length(figNames)
            figName = figNames{fIdx};
            
            % Save as .fig
            figFileName = sprintf('%s_%s.fig', baseName, figName);
            savefig(figHandles.(figName), fullfile(savePath, figFileName));
            fprintf('Saved: %s\n', figFileName);
            
            % Save as .png for easy viewing
            pngFileName = sprintf('%s_%s.png', baseName, figName);
            saveas(figHandles.(figName), fullfile(savePath, pngFileName));
            fprintf('Saved: %s\n', pngFileName);
        end
        
        fprintf('\nAll figures saved successfully!\n');
        fprintf('Total: %d figures (both .fig and .png formats)\n', length(figNames));
    end
end

disp('Analysis complete!');
