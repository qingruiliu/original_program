%% Multi-Animal Behavioral Training Analysis Script
% This script analyzes discrimination training data from multiple animals
% Load per-animal behavioral_analysis_results.mat, group comparison, and logistic fits
% Author: Modified for Liu's multi-animal behavioral analysis
% Date: 2026.05.11

clear; clc; close all;

%% 1. Get number of animals per group
prompt = {'Number of WT mice:', 'Number of APP mice:'};
dlgtitle = 'Multi-Animal Analysis';
dims = [1 40];
definput = {'3', '3'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('Analysis canceled by user');
    return;
end

nWT = str2double(answer{1});
nAPP = str2double(answer{2});

if any(isnan([nWT, nAPP])) || any([nWT, nAPP] < 0) || (nWT + nAPP) < 1
    error('Invalid number of animals');
end

groupNames = {'WT', 'APP'};
groupCounts = [nWT, nAPP];
nAnimalsTotal = sum(groupCounts);

disp(['Analyzing ' num2str(nAnimalsTotal) ' animals...']);

%% 2. Initialize storage for all animals
allAnimals = struct([]);
animalCounter = 0;

%% 3. Load and process data for each animal
for groupIdx = 1:numel(groupNames)
    genotype = groupNames{groupIdx};
    nGroup = groupCounts(groupIdx);
    
    for groupAnimalIdx = 1:nGroup
        animalCounter = animalCounter + 1;
        fprintf('\n===== Processing %s Animal %d/%d (Overall %d/%d) =====\n', ...
            genotype, groupAnimalIdx, nGroup, animalCounter, nAnimalsTotal);
        
        prompt = {sprintf('%s Animal ID:', genotype)};
        dlgtitle = sprintf('%s Animal %d Information', genotype, groupAnimalIdx);
        dims = [1 50];
        definput = {sprintf('%s%d', genotype, groupAnimalIdx)};
        animalInfo = inputdlg(prompt, dlgtitle, dims, definput);
        
        if isempty(animalInfo)
            warning('%s Animal %d skipped', genotype, groupAnimalIdx);
            continue;
        end
        
        animalID = animalInfo{1};
        
        % Select pre-processed results for this animal
        [fileName, path] = uigetfile('*.mat', ...
            sprintf('Select behavioral_analysis_results.mat for %s', animalID));
        
        if isequal(fileName, 0)
            warning('No results file selected for animal %s', animalID);
            continue;
        end
        
        fullPath = fullfile(path, fileName);
        data = load(fullPath);
        if ~isfield(data, 'metrics')
            warning('Selected file does not contain metrics: %s', fileName);
            continue;
        end
        metrics = data.metrics;
        
        if ~isfield(metrics, 'correctRate')
            warning('Missing correctRate in %s', fileName);
            continue;
        end
        
        nDays = numel(metrics.correctRate);
        fprintf('Loaded %d training days for %s...\n', nDays, animalID);
        
        % Initialize metrics for this animal
        animal = struct();
        animal.ID = animalID;
        animal.genotype = genotype;
        animal.path = path;
        animal.resultsFile = fullPath;
        if isfield(metrics, 'fileNames')
            animal.fileNames = metrics.fileNames;
        else
            animal.fileNames = {};
        end
        animal.nDays = nDays;
        
        % Metrics
        animal.correctRate = padVector(metrics.correctRate, nDays, NaN);
        animal.hitRate = padVector(getFieldOrEmpty(metrics, 'hitRate'), nDays, NaN);
        animal.FARate = padVector(getFieldOrEmpty(metrics, 'FARate'), nDays, NaN);
        animal.CRRate = padVector(getFieldOrEmpty(metrics, 'CRRate'), nDays, NaN);
        animal.missRate = padVector(getFieldOrEmpty(metrics, 'missRate'), nDays, NaN);
        
        % Trial counts
        animal.nHit = padVector(getFieldOrEmpty(metrics, 'nHit'), nDays, NaN);
        animal.nMiss = padVector(getFieldOrEmpty(metrics, 'nMiss'), nDays, NaN);
        animal.nFA = padVector(getFieldOrEmpty(metrics, 'nFA'), nDays, NaN);
        animal.nCR = padVector(getFieldOrEmpty(metrics, 'nCR'), nDays, NaN);
        animal.nTotal = padVector(getFieldOrEmpty(metrics, 'nTotal'), nDays, NaN);
        animal.nTarget = padVector(getFieldOrEmpty(metrics, 'nTarget'), nDays, NaN);
        animal.nNonTarget = padVector(getFieldOrEmpty(metrics, 'nNonTarget'), nDays, NaN);
        
        % Hit latency
        animal.hitLatency = padCell(getFieldOrEmpty(metrics, 'hitLatency'), nDays);
        animal.hitLatencyMedian = padVector(getFieldOrEmpty(metrics, 'hitLatencyMedian'), nDays, NaN);
        animal.hitLatencyStd = padVector(getFieldOrEmpty(metrics, 'hitLatencyStd'), nDays, NaN);
        
        if all(isnan(animal.missRate)) && any(~isnan(animal.hitRate))
            animal.missRate = 1 - animal.hitRate;
        end
        if all(isnan(animal.CRRate)) && any(~isnan(animal.FARate))
            animal.CRRate = 1 - animal.FARate;
        end
        if all(isnan(animal.hitLatencyMedian)) && any(~cellfun(@isempty, animal.hitLatency))
            for dayIdx = 1:nDays
                latencies = animal.hitLatency{dayIdx};
                if ~isempty(latencies)
                    animal.hitLatencyMedian(dayIdx) = median(latencies);
                    animal.hitLatencyStd(dayIdx) = std(latencies);
                end
            end
        end
        
        %% Display summary for this animal
        fprintf('\n----- Summary for %s -----\n', animalID);
        fprintf('Genotype: %s\n', genotype);
        fprintf('Total training days: %d\n', nDays);
        fprintf('Total trials: %d\n', sum(animal.nTotal, 'omitnan'));
        fprintf('Final correct rate: %.2f%%\n', animal.correctRate(end) * 100);
        fprintf('Final hit rate: %.2f%%\n', animal.hitRate(end) * 100);
        fprintf('Final FA rate: %.2f%%\n', animal.FARate(end) * 100);
        fprintf('\n');
        
        %% Create individual figure for this animal
        fig = figure('Name', sprintf('%s Analysis', animalID), 'Position', [100, 100, 1400, 900]);
        
        % Color scheme
        colorCorrect = [0.2, 0.7, 0.3];
        colorHit = [0.3, 0.5, 0.9];
        colorFA = [0.9, 0.3, 0.3];
        colorLatency = [0.5, 0.2, 0.7];
        
        days = 1:nDays;
        
        % Subplot 1: Correct Rate
        subplot(2, 3, 1);
        plot(days, animal.correctRate * 100, 'o-', 'Color', colorCorrect, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorCorrect);
        yline(80, '--k', 'LineWidth', 1.5, 'Alpha', 0.5);
        xlabel('Training Day', 'FontSize', 12);
        ylabel('Correct Rate (%)', 'FontSize', 12);
        title('Correct Rate Across Days', 'FontSize', 13);
        ylim([0, 105]);
        xlim([0.5, nDays + 0.5]);
        grid on;
        set(gca, 'FontSize', 11);
        
        % Subplot 2: Hit Rate and FA Rate
        subplot(2, 3, 2);
        hold on;
        plot(days, animal.hitRate * 100, 'o-', 'Color', colorHit, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorHit, 'DisplayName', 'Hit Rate');
        plot(days, animal.FARate * 100, 's-', 'Color', colorFA, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorFA, 'DisplayName', 'FA Rate');
        xlabel('Training Day', 'FontSize', 12);
        ylabel('Rate (%)', 'FontSize', 12);
        title('Hit Rate vs FA Rate', 'FontSize', 13);
        legend('Location', 'best', 'FontSize', 10);
        ylim([0, 105]);
        xlim([0.5, nDays + 0.5]);
        grid on;
        set(gca, 'FontSize', 11);
        hold off;
        
        % Subplot 3: Trial Counts Stacked Bar
        subplot(2, 3, 3);
        trialCounts = [animal.nHit, animal.nCR, animal.nFA, animal.nMiss];
        bar(days, trialCounts, 'stacked');
        xlabel('Training Day', 'FontSize', 12);
        ylabel('Number of Trials', 'FontSize', 12);
        title('Trial Type Distribution', 'FontSize', 13);
        legend({'Hit', 'CR', 'FA', 'Miss'}, 'Location', 'best', 'FontSize', 10);
        xlim([0.5, nDays + 0.5]);
        grid on;
        set(gca, 'FontSize', 11);
        
        % Subplot 4: Hit Latency Median with Error Bars
        subplot(2, 3, 4);
        validDays = ~isnan(animal.hitLatencyMedian);
        errorbar(days(validDays), animal.hitLatencyMedian(validDays), animal.hitLatencyStd(validDays), 'o-', ...
            'Color', colorLatency, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', colorLatency);
        xlabel('Training Day', 'FontSize', 12);
        ylabel('Hit Latency (s)', 'FontSize', 12);
        title('Hit Latency (Median ± SD)', 'FontSize', 13);
        maxLatency = max(animal.hitLatencyMedian(validDays) + animal.hitLatencyStd(validDays));
        if ~isempty(maxLatency) && maxLatency > 0
            ylim([0, maxLatency * 1.2]);
        end
        xlim([0.5, nDays + 0.5]);
        grid on;
        set(gca, 'FontSize', 11);
        
        % Subplot 5: Hit Latency Distribution
        subplot(2, 3, 5);
        hold on;
        grayColor = [0.7, 0.7, 0.7];
        for dayIdx = 1:nDays
            if ~isempty(animal.hitLatency{dayIdx})
                latencies = animal.hitLatency{dayIdx};
                x = dayIdx * ones(size(latencies));
                jitter = 0.15 * randn(size(x));
                scatter(x + jitter, latencies, 36, grayColor, 'filled', 'MarkerFaceAlpha', 0.4);
                
                % Plot median line
                medVal = animal.hitLatencyMedian(dayIdx);
                if ~isnan(medVal)
                    line([dayIdx - 0.3, dayIdx + 0.3], [medVal, medVal], 'Color', 'k', 'LineWidth', 3);
                end
            end
        end
        xlabel('Training Day', 'FontSize', 12);
        ylabel('Hit Latency (s)', 'FontSize', 12);
        title('Hit Latency Distribution', 'FontSize', 13);
        ylim([0, 5]);
        xlim([0.5, nDays + 0.5]);
        grid on;
        set(gca, 'FontSize', 11);
        hold off;
        
        % Subplot 6: Performance Matrix
        subplot(2, 3, 6);
        avgHitRate = mean(animal.hitRate, 'omitnan');
        avgMissRate = mean(animal.missRate, 'omitnan');
        avgFARate = mean(animal.FARate, 'omitnan');
        avgCRRate = mean(animal.CRRate, 'omitnan');
        
        perfMatrix = [avgHitRate, avgMissRate; avgFARate, avgCRRate] * 100;
        imagesc(perfMatrix);
        colormap(gca, flipud(gray));
        colorbar;
        set(gca, 'XTick', 1:2, 'XTickLabel', {'Lick', 'No Lick'}, ...
            'YTick', 1:2, 'YTickLabel', {'Target', 'Non-Target'}, 'FontSize', 11);
        title('Average Performance Matrix (%)', 'FontSize', 13);
        xlabel('Response', 'FontSize', 12);
        ylabel('Stimulus', 'FontSize', 12);
        
        % Add text annotations
        for i = 1:2
            for j = 1:2
                text(j, i, sprintf('%.1f%%', perfMatrix(i,j)), ...
                    'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 16, 'FontWeight', 'bold');
            end
        end
        
        sgtitle(sprintf('%s (%s) - Behavioral Training Analysis', animalID, genotype), ...
            'FontSize', 15, 'FontWeight', 'bold');
        
        % Store figure handle
        animal.figure = fig;
        
        % Store this animal's data
        if isempty(allAnimals)
            allAnimals = animal;
        else
            [allAnimals, animal] = harmonizeStructArray(allAnimals, animal);
            allAnimals(end + 1) = animal;
        end
    end
end

if isempty(allAnimals)
    error('No animals were processed.');
end

%% 4. Create comparison figure across all animals
fprintf('\n===== Creating Multi-Animal Comparison Figure =====\n');

compFig = figure('Name', 'Multi-Animal Comparison', 'Position', [150, 50, 1600, 1000]);

groupColors = struct('WT', [0.5, 0.5, 0.5], 'APP', [0.9, 0.3, 0.3]);
groupsPresent = unique({allAnimals.genotype}, 'stable');

% Subplot 1: Correct Rate Comparison (group mean ± SEM)
subplot(2, 3, 1);
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    [meanVals, semVals] = groupMeanSem(genotypeAnimals, 'correctRate');
    days = 1:numel(meanVals);
    errorbar(days, meanVals * 100, semVals * 100, 'o-', 'Color', groupColors.(genotype), ...
        'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', groupColors.(genotype), ...
        'DisplayName', genotype);
end
yline(80, '--k', 'LineWidth', 1.5, 'Alpha', 0.5);
xlabel('Training Day', 'FontSize', 12);
ylabel('Correct Rate (%)', 'FontSize', 12);
title('Correct Rate - Group Mean ± SEM', 'FontSize', 13);
ylim([0, 105]);
legend('Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 11);
hold off;

% Subplot 2: Hit Rate Comparison (group mean ± SEM)
subplot(2, 3, 2);
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    [meanVals, semVals] = groupMeanSem(genotypeAnimals, 'hitRate');
    days = 1:numel(meanVals);
    errorbar(days, meanVals * 100, semVals * 100, 'o-', 'Color', groupColors.(genotype), ...
        'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', groupColors.(genotype), ...
        'DisplayName', genotype);
end
xlabel('Training Day', 'FontSize', 12);
ylabel('Hit Rate (%)', 'FontSize', 12);
title('Hit Rate - Group Mean ± SEM', 'FontSize', 13);
ylim([0, 105]);
legend('Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 11);
hold off;

% Subplot 3: FA Rate Comparison (group mean ± SEM)
subplot(2, 3, 3);
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    [meanVals, semVals] = groupMeanSem(genotypeAnimals, 'FARate');
    days = 1:numel(meanVals);
    errorbar(days, meanVals * 100, semVals * 100, 's-', 'Color', groupColors.(genotype), ...
        'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', groupColors.(genotype), ...
        'DisplayName', genotype);
end
xlabel('Training Day', 'FontSize', 12);
ylabel('FA Rate (%)', 'FontSize', 12);
title('False Alarm Rate - Group Mean ± SEM', 'FontSize', 13);
ylim([0, 105]);
legend('Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 11);
hold off;

% Subplot 4: Hit Latency Comparison (group mean ± SEM)
subplot(2, 3, 4);
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    [meanVals, semVals] = groupMeanSem(genotypeAnimals, 'hitLatencyMedian');
    days = 1:numel(meanVals);
    errorbar(days, meanVals, semVals, 'o-', 'Color', groupColors.(genotype), ...
        'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', groupColors.(genotype), ...
        'DisplayName', genotype);
end
xlabel('Training Day', 'FontSize', 12);
ylabel('Hit Latency (s)', 'FontSize', 12);
title('Hit Latency - Group Mean ± SEM', 'FontSize', 13);
ylim([0, 4]);
legend('Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 11);
hold off;

% Subplot 5: Average Performance by Genotype
subplot(2, 3, 5);
nGenotypes = numel(groupsPresent);
avgPerf = zeros(nGenotypes, 3); % Correct, Hit, FA rates
semPerf = zeros(nGenotypes, 3);

for gIdx = 1:nGenotypes
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    
    allCorrect = [];
    allHit = [];
    allFA = [];
    
    for aIdx = 1:length(genotypeAnimals)
        allCorrect = [allCorrect; genotypeAnimals(aIdx).correctRate];
        allHit = [allHit; genotypeAnimals(aIdx).hitRate];
        allFA = [allFA; genotypeAnimals(aIdx).FARate];
    end
    
    avgPerf(gIdx, 1) = mean(allCorrect, 'omitnan') * 100;
    avgPerf(gIdx, 2) = mean(allHit, 'omitnan') * 100;
    avgPerf(gIdx, 3) = mean(allFA, 'omitnan') * 100;
    
    semPerf(gIdx, 1) = std(allCorrect, 0, 'omitnan') / sqrt(sum(~isnan(allCorrect))) * 100;
    semPerf(gIdx, 2) = std(allHit, 0, 'omitnan') / sqrt(sum(~isnan(allHit))) * 100;
    semPerf(gIdx, 3) = std(allFA, 0, 'omitnan') / sqrt(sum(~isnan(allFA))) * 100;
end

bar(avgPerf);
set(gca, 'XTickLabel', groupsPresent, 'FontSize', 11);
xlabel('Genotype', 'FontSize', 12);
ylabel('Rate (%)', 'FontSize', 12);
title('Average Performance by Genotype', 'FontSize', 13);
legend({'Correct Rate', 'Hit Rate', 'FA Rate'}, 'Location', 'best', 'FontSize', 10);
ylim([0, 105]);
grid on;

% Subplot 6: Final Day Performance Comparison
subplot(2, 3, 6);
finalPerf = zeros(nGenotypes, 3);
finalSem = zeros(nGenotypes, 3);

for gIdx = 1:nGenotypes
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    finalCorrect = arrayfun(@(a) a.correctRate(end) * 100, genotypeAnimals);
    finalHit = arrayfun(@(a) a.hitRate(end) * 100, genotypeAnimals);
    finalFA = arrayfun(@(a) a.FARate(end) * 100, genotypeAnimals);
    
    finalPerf(gIdx, :) = [mean(finalCorrect, 'omitnan'), mean(finalHit, 'omitnan'), mean(finalFA, 'omitnan')];
    finalSem(gIdx, :) = [std(finalCorrect, 0, 'omitnan'), std(finalHit, 0, 'omitnan'), std(finalFA, 0, 'omitnan')] ...
        ./ sqrt([sum(~isnan(finalCorrect)), sum(~isnan(finalHit)), sum(~isnan(finalFA))]);
end

bar(finalPerf);
set(gca, 'XTickLabel', groupsPresent, 'FontSize', 11);
xlabel('Genotype', 'FontSize', 12);
ylabel('Rate (%)', 'FontSize', 12);
title('Final Day Performance by Genotype', 'FontSize', 13);
legend({'Correct', 'Hit', 'FA'}, 'Location', 'best', 'FontSize', 10);
ylim([0, 105]);
grid on;

sgtitle('Multi-Animal Behavioral Analysis Comparison', 'FontSize', 16, 'FontWeight', 'bold');

%% 5. Logistic regression fitting (learning curves)
fprintf('\n===== Logistic Regression Fitting =====\n');

hasLsqcurvefit = exist('lsqcurvefit', 'file') == 2;
if ~hasLsqcurvefit
    warning('lsqcurvefit not found. Logistic fitting skipped.');
else
    for idx = 1:length(allAnimals)
        [fitParams, fitDays, fitCurve] = fitLogisticCurve(allAnimals(idx).correctRate);
        allAnimals(idx).logisticParams = fitParams;
        allAnimals(idx).logisticFitDays = fitDays;
        allAnimals(idx).logisticFitCurve = fitCurve;
    end
    
    % Extract fit parameters for group comparison
    slopeVals = struct();
    inflectionVals = struct();
    upperVals = struct();
    lowerVals = struct();
    
    for gIdx = 1:nGenotypes
        genotype = groupsPresent{gIdx};
        genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
        params = vertcat(genotypeAnimals.logisticParams);
        params = params(~any(isnan(params), 2), :);
        
        slopeVals.(genotype) = params(:, 3);
        inflectionVals.(genotype) = params(:, 4);
        lowerVals.(genotype) = params(:, 1);
        upperVals.(genotype) = params(:, 2);
    end
    
    % Comparison figure
    fitFig = figure('Name', 'Learning Curve Logistic Fits', 'Position', [200, 80, 1500, 900]);
    
    subplot(2, 2, 1);
    hold on;
    for idx = 1:length(allAnimals)
        animal = allAnimals(idx);
        days = 1:animal.nDays;
        color = groupColors.(animal.genotype);
        plot(days, animal.correctRate * 100, 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 6);
        if ~isempty(animal.logisticFitDays)
            plot(animal.logisticFitDays, animal.logisticFitCurve * 100, '-', 'Color', color, 'LineWidth', 2);
        end
    end
    xlabel('Training Day', 'FontSize', 12);
    ylabel('Correct Rate (%)', 'FontSize', 12);
    title('Per-Animal Logistic Fits', 'FontSize', 13);
    ylim([0, 105]);
    grid on;
    set(gca, 'FontSize', 11);
    hold off;
    
    subplot(2, 2, 2);
    boxplotData = [];
    boxplotGroup = {};
    for gIdx = 1:nGenotypes
        genotype = groupsPresent{gIdx};
        vals = slopeVals.(genotype);
        boxplotData = [boxplotData; vals];
        boxplotGroup = [boxplotGroup; repmat({genotype}, numel(vals), 1)];
    end
    boxplot(boxplotData, boxplotGroup);
    ylabel('Slope (k)', 'FontSize', 12);
    title('Learning Rate (Slope)', 'FontSize', 13);
    grid on;
    
    subplot(2, 2, 3);
    boxplotData = [];
    boxplotGroup = {};
    for gIdx = 1:nGenotypes
        genotype = groupsPresent{gIdx};
        vals = inflectionVals.(genotype);
        boxplotData = [boxplotData; vals];
        boxplotGroup = [boxplotGroup; repmat({genotype}, numel(vals), 1)];
    end
    boxplot(boxplotData, boxplotGroup);
    ylabel('Inflection Day (x0)', 'FontSize', 12);
    title('Learning Onset (Inflection)', 'FontSize', 13);
    grid on;
    
    subplot(2, 2, 4);
    boxplotData = [];
    boxplotGroup = {};
    for gIdx = 1:nGenotypes
        genotype = groupsPresent{gIdx};
        vals = upperVals.(genotype);
        boxplotData = [boxplotData; vals];
        boxplotGroup = [boxplotGroup; repmat({genotype}, numel(vals), 1)];
    end
    boxplot(boxplotData, boxplotGroup);
    ylabel('Upper Asymptote', 'FontSize', 12);
    title('Final Performance (Upper)', 'FontSize', 13);
    grid on;
    
    % Statistical comparison (WT vs APP)
    if isfield(slopeVals, 'WT') && isfield(slopeVals, 'APP')
        if numel(slopeVals.WT) >= 2 && numel(slopeVals.APP) >= 2
            [~, pSlope] = ttest2(slopeVals.WT, slopeVals.APP);
            [~, pX0] = ttest2(inflectionVals.WT, inflectionVals.APP);
            fprintf('Slope comparison WT vs APP: p = %.4f\n', pSlope);
            fprintf('Inflection comparison WT vs APP: p = %.4f\n', pX0);
        else
            fprintf('Not enough animals for statistical comparison.\n');
        end
    end
end

%% 6. Save all results
answer = questdlg('Do you want to save all analysis results?', 'Save Results', 'Yes', 'No', 'Yes');

if strcmp(answer, 'Yes')
    [saveFileName, savePath] = uiputfile('*.mat', 'Save multi-animal analysis results', 'multi_animal_results.mat');
    
    if ~isequal(saveFileName, 0)
        % Save data
        save(fullfile(savePath, saveFileName), 'allAnimals');
        fprintf('Results saved to: %s\n', fullfile(savePath, saveFileName));
        
        % Save individual animal figures
        for idx = 1:length(allAnimals)
            animalFigName = sprintf('%s_analysis.fig', allAnimals(idx).ID);
            savefig(allAnimals(idx).figure, fullfile(savePath, animalFigName));
            fprintf('Figure saved: %s\n', animalFigName);
        end
        
        % Save comparison figure
        compFigName = 'multi_animal_comparison.fig';
        savefig(compFig, fullfile(savePath, compFigName));
        fprintf('Comparison figure saved: %s\n', compFigName);
    end
end

fprintf('\n===== Multi-Animal Analysis Complete! =====\n');
fprintf('Analyzed %d animals\n', length(allAnimals));
fprintf('Individual figures and comparison figure created\n');

%% Local functions
function [meanVals, semVals] = groupMeanSem(groupAnimals, fieldName)
    maxDays = max([groupAnimals.nDays]);
    data = NaN(length(groupAnimals), maxDays);
    for idx = 1:length(groupAnimals)
        vals = groupAnimals(idx).(fieldName);
        data(idx, 1:length(vals)) = vals(:)';
    end
    meanVals = mean(data, 1, 'omitnan');
    semVals = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
end

function [fitParams, fitDays, fitCurve] = fitLogisticCurve(rateVector)
    y = rateVector(:);
    x = (1:numel(y))';
    valid = ~isnan(y);
    x = x(valid);
    y = y(valid);
    
    if numel(x) < 3
        fitParams = nan(1, 4);
        fitDays = [];
        fitCurve = [];
        return;
    end
    
    lb = [0, 0, 0, 1];
    ub = [1, 1, 5, max(x)];
    p0 = [max(0, min(y)), min(1, max(y)), 1, median(x)];
    
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    fitParams = lsqcurvefit(@logisticFun, p0, x, y, lb, ub, options);
    
    fitDays = linspace(min(x), max(x), 100)';
    fitCurve = logisticFun(fitParams, fitDays);
end

function y = logisticFun(p, x)
    lower = p(1);
    upper = p(2);
    slope = p(3);
    x0 = p(4);
    y = lower + (upper - lower) ./ (1 + exp(-slope * (x - x0)));
end

function vals = getFieldOrEmpty(s, fieldName)
    if isstruct(s) && isfield(s, fieldName)
        vals = s.(fieldName);
    else
        vals = [];
    end
end

function out = padVector(vals, nDays, fillValue)
    if nargin < 3
        fillValue = NaN;
    end
    out = fillValue * ones(nDays, 1);
    if isempty(vals)
        return;
    end
    vals = vals(:);
    n = min(nDays, numel(vals));
    out(1:n) = vals(1:n);
end

function out = padCell(vals, nDays)
    out = cell(nDays, 1);
    if isempty(vals)
        return;
    end
    n = min(nDays, numel(vals));
    out(1:n) = vals(1:n);
end

function [structArray, newStruct] = harmonizeStructArray(structArray, newStruct)
    arrayFields = fieldnames(structArray);
    newFields = fieldnames(newStruct);
    missingInArray = setdiff(newFields, arrayFields);
    for idx = 1:numel(missingInArray)
        [structArray.(missingInArray{idx})] = deal([]);
    end
    missingInNew = setdiff(arrayFields, newFields);
    for idx = 1:numel(missingInNew)
        newStruct.(missingInNew{idx}) = [];
    end
end