%% Multi-animal summary figures from multi_animal_results.mat
% Creates and saves 4 standalone figures from allAnimals.

clear; clc; close all;

%% 1. Load results
[fileName, path] = uigetfile('*.mat', 'Select multi_animal_results.mat');
if isequal(fileName, 0)
    disp('User canceled file selection');
    return;
end

fullPath = fullfile(path, fileName);
data = load(fullPath);
if ~isfield(data, 'allAnimals')
    error('Selected file does not contain allAnimals');
end

allAnimals = data.allAnimals;
if isempty(allAnimals)
    error('allAnimals is empty');
end

%% 2. Common settings
colors = struct('WT', [0, 0, 0], 'APP', [1, 0.5, 0]);
groupsPresent = unique({allAnimals.genotype}, 'stable');

outDir = path;

%% Figure 1: Correct rate comparison (line)
fig1 = figure('Name', 'Correct Rate Comparison', 'Position', [100, 100, 700, 550]);
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    % Individual animal traces (light color)
    lineColor = colors.(genotype) * 0.35 + 0.65;
    for aIdx = 1:numel(genotypeAnimals)
        days = 1:genotypeAnimals(aIdx).nDays;
        plot(days, genotypeAnimals(aIdx).correctRate * 100, '-', ...
            'Color', lineColor, 'LineWidth', 4, 'HandleVisibility', 'off');
    end
end

yline(80, '--', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 2, 'Alpha', 0.5);
xlabel('Training Day', 'FontSize', 16);
ylabel('Correct Rate (%)', 'FontSize', 16);
title('Correct Rate - Group Mean ± SEM', 'FontSize', 16);
title('Correct Rate per Animal', 'FontSize', 16);
ylim([40, 100]);
grid off;
set(gca, 'FontSize', 18);
hold off;

savefig(fig1, fullfile(outDir, 'fig1_correct_rate_comparison.fig'));

%% Figure 2: Sessions to expert (>=80% correct)
fig2 = figure('Name', 'Sessions to Expert', 'Position', [120, 120, 600, 520]);
APP = [9,15,9,10,7,5];
WT = [10,9,7,13,9];
groupData = struct('WT', WT(:), 'APP', APP(:));
meanSessions = zeros(numel(groupsPresent), 1);
semSessions = zeros(numel(groupsPresent), 1);
barColors = zeros(numel(groupsPresent), 3);
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    daysToExpert = groupData.(genotype);
    meanSessions(gIdx) = mean(daysToExpert, 'omitnan');
    semSessions(gIdx) = std(daysToExpert, 0, 'omitnan') / sqrt(sum(~isnan(daysToExpert)));
    if strcmp(genotype, 'WT')
        barColors(gIdx, :) = [0.6, 0.6, 0.6];
    else
        barColors(gIdx, :) = colors.(genotype);
    end
end

bar(meanSessions, 'FaceColor', 'flat');
set(get(gca, 'Children'), 'CData', barColors);

hold on;
errorbar(1:numel(groupsPresent), meanSessions, semSessions, 'k.', 'LineWidth', 1.5, 'CapSize', 12);

% WT vs APP p-value (t-test)
pText = 'p = n/a';
if numel(WT) >= 2 && numel(APP) >= 2
    [~, pVal] = ttest2(WT, APP);
    pText = sprintf('p = %.3g', pVal);
end
yMax = max(meanSessions + semSessions);
text(1.5, yMax * 1.08, pText, 'HorizontalAlignment', 'center', 'FontSize', 12);

set(gca, 'XTick', 1:numel(groupsPresent), 'XTickLabel', groupsPresent, 'FontSize', 11);
xlabel('Genotype', 'FontSize', 12);
ylabel('Sessions to Expert', 'FontSize', 16);
title('Sessions to Reach >=80% Correct', 'FontSize', 13);
grid off;
box off;
hold off;
set(gca, 'FontSize', 18);
savefig(fig2, fullfile(outDir, 'fig2_sessions_to_expert.fig'));

%% Figure 3: Quartile day hit latency (mean ± SEM)
fig3 = figure('Name', 'Quartile Hit Latency', 'Position', [140, 140, 700, 550]);
hold on;

quartilePoints = [0.25, 0.5, 0.75, 1.0];
quartileLabels = {'25%', '50%', '75%', '100%'};

for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    [meanVals, semVals] = groupQuartileLatency(genotypeAnimals, quartilePoints);
    errorbar(1:numel(quartilePoints), meanVals, semVals, 'o-', ...
        'Color', colors.(genotype), 'LineWidth', 2.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', colors.(genotype), 'DisplayName', genotype);
end

set(gca, 'XTick', 1:numel(quartilePoints), 'XTickLabel', quartileLabels, 'FontSize', 11);
xlabel('Training Progress', 'FontSize', 12);
ylabel('Hit Trial Lick Latency (s)', 'FontSize', 12);
title('Hit Latency at Quartile Days (Mean ± SEM)', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
grid off;
hold off;

set(gca, 'FontSize', 18);

savefig(fig3, fullfile(outDir, 'fig3_quartile_hit_latency.fig'));

%% Figure 4a: Per-animal d' scatter + logistic fits
fig4a = figure('Name', 'Per-Animal Fits', 'Position', [160, 160, 750, 550]);
hold on;
maxDays = max([allAnimals.nDays]);
fitDays = (1:maxDays)';

for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    lineColor = colors.(genotype) * 0.35 + 0.65;
    for aIdx = 1:numel(genotypeAnimals)
        days = 1:genotypeAnimals(aIdx).nDays;
        dprimeVals = getDprimeVector(genotypeAnimals(aIdx));
        scatter(days, dprimeVals, 18, ...
            'MarkerFaceColor', lineColor, 'MarkerEdgeColor', lineColor, 'MarkerFaceAlpha', 0.6);
        [params, valid] = getDprimeLogisticParams(genotypeAnimals(aIdx));
        if valid
            curve = logisticFun(params, fitDays);
            plot(fitDays, curve, '-', 'Color', lineColor, 'LineWidth', 1.2);
        end
    end
end

xlabel('Training Day', 'FontSize', 12);
ylabel('d''', 'FontSize', 12);
title('Per-Animal d'' + Logistic Fits', 'FontSize', 13);
grid on;
hold off;

savefig(fig4a, fullfile(outDir, 'fig4a_per_animal_scatter_fits.fig'));

%% Figure 4b: Mean ± SEM d' logistic fits + slope comparison
figure;
hold on;
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    curves = [];
    for aIdx = 1:numel(genotypeAnimals)
        [params, valid] = getDprimeLogisticParams(genotypeAnimals(aIdx));
        if valid
            curve = logisticFun(params, fitDays);
            curves = [curves, curve];
        end
    end
    if ~isempty(curves)
        meanCurve = mean(curves, 2, 'omitnan');
        semCurve = std(curves, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(curves), 2));
        fill([fitDays; flipud(fitDays)], ...
            [(meanCurve - semCurve); flipud((meanCurve + semCurve))], ...
            colors.(genotype), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        plot(fitDays, meanCurve, '-', 'Color', colors.(genotype), 'LineWidth', 2.5, ...
            'DisplayName', genotype);
    end
end

xlabel('Training Day', 'FontSize', 12);
ylabel('Discriminability Index (d'')', 'FontSize', 12);
title('Group d'' Fits (Mean ± SEM)', 'FontSize', 13);
legend off;
ylim([-0.5 6]);
set(gca, 'FontSize', 18);
grid off;
hold off;

figure;
slopes = cell(numel(groupsPresent), 1);
for gIdx = 1:numel(groupsPresent)
    genotype = groupsPresent{gIdx};
    genotypeAnimals = allAnimals(strcmp({allAnimals.genotype}, genotype));
    slopeVals = [];
    for aIdx = 1:numel(genotypeAnimals)
        [params, valid] = getDprimeLogisticParams(genotypeAnimals(aIdx));
        if valid
            slopeVals = [slopeVals; params(3)];
        end
    end
    slopes{gIdx} = slopeVals;
end

meanSlopes = zeros(numel(groupsPresent), 1);
semSlopes = zeros(numel(groupsPresent), 1);
barColors = zeros(numel(groupsPresent), 3);
for gIdx = 1:numel(groupsPresent)
    vals = slopes{gIdx};
    meanSlopes(gIdx) = mean(vals, 'omitnan');
    semSlopes(gIdx) = std(vals, 0, 'omitnan') / sqrt(sum(~isnan(vals)));
    if strcmp(groupsPresent{gIdx}, 'WT')
        barColors(gIdx, :) = [0.6, 0.6, 0.6];
    else
        barColors(gIdx, :) = colors.(groupsPresent{gIdx});
    end
end

bar(meanSlopes, 'FaceColor', 'flat');
set(get(gca, 'Children'), 'CData', barColors);
hold on;
errorbar(1:numel(groupsPresent), meanSlopes, semSlopes, 'k.', 'LineWidth', 1.5, 'CapSize', 12);
for gIdx = 1:numel(groupsPresent)
    vals = slopes{gIdx};
    if isempty(vals)
        continue;
    end
    xJitter = gIdx + 0.12 * randn(size(vals));
    scatter(xJitter, vals, 30, colors.(groupsPresent{gIdx}), 'filled', 'MarkerFaceAlpha', 0.6);
end

% WT vs APP p-value (t-test)
pText = 'p = n/a';
if numel(slopes) == 2 && numel(slopes{1}) >= 2 && numel(slopes{2}) >= 2
    [~, pVal] = ttest2(slopes{1}, slopes{2});
    pText = sprintf('p = %.3g', pVal);
end
yMax = max(cellfun(@max, slopes));
text(1.5, yMax * 1.08, pText, 'HorizontalAlignment', 'center', 'FontSize', 11);

ylabel('Slope (k)', 'FontSize', 12);
title('Learning Rate (Slope)', 'FontSize', 13);
grid off;
hold off;
box off;
set(gca, 'FontSize', 18);


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

function day = firstExpertDay(rateVector, threshold)
    idx = find(rateVector >= threshold, 1, 'first');
    if isempty(idx)
        day = NaN;
    else
        day = idx;
    end
end

function [meanVals, semVals] = groupQuartileLatency(groupAnimals, quartilePoints)
    nQ = numel(quartilePoints);
    vals = NaN(numel(groupAnimals), nQ);
    for aIdx = 1:numel(groupAnimals)
        animal = groupAnimals(aIdx);
        nDays = animal.nDays;
        qDays = unique(max(1, min(nDays, ceil(quartilePoints * nDays))));
        qDays = padArrayToLength(qDays, nQ, nDays);
        for qIdx = 1:nQ
            dayIdx = qDays(qIdx);
            latency = NaN;
            if isfield(animal, 'hitLatencyMedian') && numel(animal.hitLatencyMedian) >= dayIdx
                latency = animal.hitLatencyMedian(dayIdx);
            elseif isfield(animal, 'hitLatency') && numel(animal.hitLatency) >= dayIdx
                latencies = animal.hitLatency{dayIdx};
                if ~isempty(latencies)
                    latency = median(latencies);
                end
            end
            vals(aIdx, qIdx) = latency;
        end
    end
    meanVals = mean(vals, 1, 'omitnan');
    semVals = std(vals, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(vals), 1));
end

function out = padArrayToLength(vals, targetLen, fillValue)
    out = fillValue * ones(1, targetLen);
    n = min(targetLen, numel(vals));
    out(1:n) = vals(1:n);
end

function dprimeVals = getDprimeVector(animal)
    if ~isfield(animal, 'hitRate') || ~isfield(animal, 'FARate')
        dprimeVals = nan(animal.nDays, 1);
        return;
    end
    hitRate = animal.hitRate(:);
    faRate = animal.FARate(:);
    epsVal = 1e-4;
    hitRate = min(max(hitRate, epsVal), 1 - epsVal);
    faRate = min(max(faRate, epsVal), 1 - epsVal);
    dprimeVals = norminv(hitRate) - norminv(faRate);
end

function [params, valid] = getDprimeLogisticParams(animal)
    valid = false;
    params = nan(1, 4);
    dprimeVals = getDprimeVector(animal);
    if all(isnan(dprimeVals))
        return;
    end
    [params, ~, ~] = fitLogisticCurveDprime(dprimeVals);
    valid = ~any(isnan(params));
end

function [fitParams, fitDays, fitCurve] = fitLogisticCurveDprime(dprimeVals)
    y = dprimeVals(:);
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
    
    yMin = min(y);
    yMax = max(y);
    if yMin == yMax
        fitParams = nan(1, 4);
        fitDays = [];
        fitCurve = [];
        return;
    end
    
    lb = [yMin - 1, yMin, 0, 1];
    ub = [yMax, yMax + 1, 5, max(x)];
    p0 = [yMin, yMax, 1, median(x)];
    
    options = optimoptions('lsqcurvefit', 'Display', 'off');
    fitParams = lsqcurvefit(@logisticFun, p0, x, y, lb, ub, options);
    
    fitDays = linspace(min(x), max(x), 100)';
    fitCurve = logisticFun(fitParams, fitDays);
end

function [params, valid] = getLogisticParams(animal)
    valid = false;
    params = nan(1, 4);
    if isfield(animal, 'logisticParams') && numel(animal.logisticParams) == 4
        params = animal.logisticParams;
        valid = ~any(isnan(params));
        return;
    end
    if isfield(animal, 'correctRate')
        [params, ~, ~] = fitLogisticCurve(animal.correctRate);
        valid = ~any(isnan(params));
    end
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

function labels = groupLabelsFromCell(dataCell, groupNames)
    labels = {};
    for idx = 1:numel(dataCell)
        labels = [labels; repmat(groupNames(idx), numel(dataCell{idx}), 1)];
    end
end
