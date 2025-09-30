%% Part 1: choose the demo .mat file
[fileName, path] = uigetfile('*.mat', 'Select the demo data for analysis', 'MultiSelect', 'off');
if isequal(fileName, 0)
    disp('User canceled the file selection.');
    return;
end
cd(path);
% Load the selected .mat file
load(fileName);

% save the actual center and the transformed center of each ROI
center = ROISegTraceTable.Center;
transCenter = ROISegTraceTable.TransformedCenter;
idx = ROISegTraceTable.ROIIndex;

% 计算每对ROI的距离
pair_idx1 = [];
pair_idx2 = [];
pair_dist = [];
pair_tanDist = [];
for i = 1:size(transCenter,1)
    for j = i+1:size(transCenter,1)
        dist = norm(center(i,:) - center(j,:));
        tanDist = norm(transCenter(i,1:2) - transCenter(j,1:2));
        if tanDist <= 10 && dist > 20
            pair_idx1 = [pair_idx1; idx(i)];
            pair_idx2 = [pair_idx2; idx(j)];
            pair_dist = [pair_dist; dist];
            pair_tanDist = [pair_tanDist; tanDist];
        end
    end
end
% Create a table to store the aligned pairs and their distances
aligned_pairs_table = table(pair_idx1, pair_idx2, pair_dist, pair_tanDist, ...
    'VariableNames', {'ROI1', 'ROI2', 'Dist', 'TanDist'});
aligned_pairs_table.R = cell(height(aligned_pairs_table), 1);
aligned_pairs_table.Features = cell(height(aligned_pairs_table), 1);

% 针对每一对aligned_pairs_table中的ROI对，进行分析并保存结果
num_pairs = height(aligned_pairs_table);
fs = 10;
dt = 1/fs;
low_cutoff = 0.1;
high_cutoff = 1.15;
[b, a] = butter(2, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

wb = waitbar(0, 'Processing features for each pair of ROIs...');
% 遍历每一对ROI
for p = 1:num_pairs
    waitbar(p/num_pairs, wb, sprintf('Processing pair %d of %d', p, num_pairs));
    roi1 = aligned_pairs_table.ROI1(p);
    roi2 = aligned_pairs_table.ROI2(p);
    cell1 = ROISegTraceTable.Segmented_trace{find(ROISegTraceTable.ROIIndex == roi1,1)};
    cell2 = ROISegTraceTable.Segmented_trace{find(ROISegTraceTable.ROIIndex == roi2,1)};
    
    % --- START: Manual Signal Contamination Check ---
    cell1_totalTrace = vertcat(cell1.traces{:}); 
    cell1_totalTrace = cell1_totalTrace(:, 3); % 取Z-score
    cell2_totalTrace = vertcat(cell2.traces{:});
    cell2_totalTrace = cell2_totalTrace(:, 3); % 取Z-score

    min_length = min(length(cell1_totalTrace), length(cell2_totalTrace));
    cell1_totalTrace = cell1_totalTrace(1:min_length);
    cell2_totalTrace = cell2_totalTrace(1:min_length);
    tempR = corrcoef(cell1_totalTrace, cell2_totalTrace);
    tempR = tempR(1,2);

    figure('Name','Signal Contamination Check');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'FontSize', 20);
    scatter(cell1_totalTrace, cell2_totalTrace, 36, 'filled', 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor', 'none');
    xlabel(sprintf('ROI%d Z-score', roi1), 'FontSize', 20);
    ylabel(sprintf('ROI%d Z-score', roi2), 'FontSize', 20);
    title(sprintf('Scatter Plot of ROI%d and ROI%d, R = %.2f', roi1, roi2, tempR), 'FontSize', 20);

    out = questdlg('Is there signal contamination in this pair?', ...
        'Signal Contamination Check', 'Yes', 'No', 'No');
    if strcmp(out, 'Yes')
        aligned_pairs_table.R{p} = NaN;
        aligned_pairs_table.Features{p} = []; % 标记为空，以便后续跳过
        close(gcf);
        continue; % 跳过后续处理
    elseif strcmp(out, 'No')
        aligned_pairs_table.R{p} = tempR;
        close(gcf); % 关闭检查窗口，继续执行
    else % User closed the dialog
        close(gcf);
        continue;
    end
    % --- END: Manual Signal Contamination Check ---

    trialNum = height(cell1.traces);
    cell1_trace = cell1.traces(1:trialNum);
    cell2_trace = cell2.traces(1:trialNum);
    
    cell1_interp = cell(trialNum, 1);
    cell2_interp = cell(trialNum, 1);
    
    for i = 1:trialNum
        min_time = min([cell1_trace{i}(:, 2); cell2_trace{i}(:, 2)]);
        max_time = max([cell1_trace{i}(:, 2); cell2_trace{i}(:, 2)]);
        common_axis = min_time:dt:max_time;
        cell1_interp{i} = [common_axis', interp1(cell1_trace{i}(:, 2), cell1_trace{i}(:, 3), common_axis, 'linear', 'extrap')'];
        cell2_interp{i} = [common_axis', interp1(cell2_trace{i}(:, 2), cell2_trace{i}(:, 3), common_axis, 'linear', 'extrap')'];
    end
    
    max_len = max(cellfun(@(x) size(x,1), cell1_interp));
    aligned_traces1 = NaN(max_len, trialNum);
    aligned_traces2 = NaN(max_len, trialNum);
    for i = 1:trialNum
        len = size(cell1_interp{i}, 1);
        aligned_traces1(1:len, i) = cell1_interp{i}(:, 2);
        aligned_traces2(1:len, i) = cell2_interp{i}(:, 2);
    end
    mean_trace1 = nanmean(aligned_traces1, 2);
    mean_trace2 = nanmean(aligned_traces2, 2);

    PLV = zeros(trialNum, 1);
    SignalCorr = zeros(trialNum, 1);
    NoiseCorr = zeros(trialNum, 1);
    EventCoincidence = zeros(trialNum, 1);
    MutualInformation = zeros(trialNum, 1);

    for i = 1:trialNum
        analysis_window_idx = cell1_interp{i}(:, 1) >= 0 & cell1_interp{i}(:, 1) <= 5;
        trace1_win = cell1_interp{i}(analysis_window_idx, 2);
        trace2_win = cell2_interp{i}(analysis_window_idx, 2);
        
        if length(trace1_win) < 2 || any(isnan(trace1_win)) || any(isnan(trace2_win))
            PLV(i)=NaN; SignalCorr(i)=NaN; NoiseCorr(i)=NaN; EventCoincidence(i)=NaN; MutualInformation(i)=NaN;
            continue;
        end

        SignalCorr(i) = corr(trace1_win, trace2_win);

        len_win = length(trace1_win);
        mean1_win = mean_trace1(analysis_window_idx);
        mean2_win = mean_trace2(analysis_window_idx);
        if length(mean1_win) == len_win && length(mean2_win) == len_win
            residual1 = trace1_win - mean1_win;
            residual2 = trace2_win - mean2_win;
            NoiseCorr(i) = corr(residual1, residual2);
        else
            NoiseCorr(i) = NaN;
        end

        event_thresh = 2.0;
        events1 = trace1_win > event_thresh;
        events2 = trace2_win > event_thresh;
        coincident_events = sum(events1 & events2);
        total_possible_events = sum(events1 | events2);
        if total_possible_events > 0
            EventCoincidence(i) = coincident_events / total_possible_events;
        else
            EventCoincidence(i) = 0;
        end

        num_bins = 10;
        combined_traces = [trace1_win; trace2_win];
        [~, edges] = histcounts(combined_traces, num_bins);
        binned1 = discretize(trace1_win, edges);
        binned2 = discretize(trace2_win, edges);
        if any(isnan(binned1)) || any(isnan(binned2))
            MutualInformation(i) = NaN;
        else
            MutualInformation(i) = mi(binned1, binned2);
        end

        cell1_filtered = filtfilt(b, a, trace1_win);
        cell2_filtered = filtfilt(b, a, trace2_win);
        phase_diff = angle(hilbert(cell1_filtered)) - angle(hilbert(cell2_filtered));
        PLV(i) = abs(mean(exp(1i * phase_diff)));
    end
    
    result = cell1.trialResult(1:trialNum);
    contrast = cell1.trialContrast(1:trialNum);
    
    feature_data = table(result, contrast, PLV, SignalCorr, NoiseCorr, EventCoincidence, MutualInformation, ...
        'VariableNames', {'Result', 'Contrast', 'PLV', 'SignalCorr', 'NoiseCorr', 'EventCoincidence', 'MutualInformation'});
    aligned_pairs_table.Features{p} = feature_data;
end
close(wb);
disp('Feature extraction completed.');
%% ==================================================
%  Part 2: Decoding Analysis
%  ==================================================
disp('Starting decoding analysis...');

% 1. 准备数据集
all_features_matrix = [];
all_behavioral_labels = [];

for p = 1:height(aligned_pairs_table)
    if isempty(aligned_pairs_table.Features{p})
        continue;
    end
    pair_features_table = aligned_pairs_table.Features{p};
    pair_features_table = rmmissing(pair_features_table);
    if isempty(pair_features_table)
        continue;
    end
    
    features_for_this_pair = [pair_features_table.PLV, pair_features_table.SignalCorr, ...
                              pair_features_table.NoiseCorr, pair_features_table.EventCoincidence, ...
                              pair_features_table.MutualInformation];
    
    decision_labels = pair_features_table.Result;
    contrast_labels = pair_features_table.Contrast;

    % << 核心修改：根据您的规则创建新的三状态行为标签 >>
    % 初始化标签向量
    num_trials_in_pair = length(decision_labels);
    behavioral_state_labels = zeros(num_trials_in_pair, 1);

    % 条件1: 正确感知与行动 (Correct Perception & Action)
    idx_state1 = (contrast_labels >= 0.1) & (decision_labels == 1 | decision_labels == 4);
    behavioral_state_labels(idx_state1) = 1;

    % 条件2: 感知到但行动错误 (Perceived but Incorrect Action)
    idx_state2 = (contrast_labels >= 0.1) & (decision_labels == 2 | decision_labels == 3);
    behavioral_state_labels(idx_state2) = 2;

    % 条件3: 猜测/未感知 (Guessing / No Perception)
    idx_state3 = contrast_labels < 0.1;
    behavioral_state_labels(idx_state3) = 3;

    % 汇总所有数据
    all_features_matrix = [all_features_matrix; features_for_this_pair];
    all_behavioral_labels = [all_behavioral_labels; behavioral_state_labels];
end

if isempty(all_features_matrix)
    disp('No valid data available for decoding. Exiting.');
    return;
end

% 2. 训练和评估分类器
disp('Training and evaluating models for different features...');

% 定义要测试的特征集
feature_sets = { ...
    {'PLV'},                all_features_matrix(:, 1); ...
    {'SignalCorr'},         all_features_matrix(:, 2); ...
    {'NoiseCorr'},          all_features_matrix(:, 3); ...
    {'EventCoincidence'},   all_features_matrix(:, 4); ...
    {'MutualInformation'},  all_features_matrix(:, 5); ...
    {'All Features'},       all_features_matrix(:, 1:5) ...
};

Y = all_behavioral_labels; % 三种状态的标签

if length(unique(Y)) < 2
    disp('Not enough classes in the data to perform classification.');
    return;
end

fprintf('============================================================\n');
fprintf('Decoding Result (Correct vs. Incorrect Choice)\n');
fprintf('Model: Logistic Regression with 10-fold Cross-Validation\n');
fprintf('------------------------------------------------------------\n');

% 循环遍历每个特征集
for f = 1:size(feature_sets, 1)
    feature_name = feature_sets{f, 1};
    X = feature_sets{f, 2};
    
    try
        cv = cvpartition(Y, 'KFold', 10);
    catch ME
        fprintf('Could not test feature set "%s". Error: %s\n', feature_name{1}, ME.message);
        continue;
    end

    accuracy_sum = 0;
    confMat_total = zeros(numel(unique(Y)), numel(unique(Y)));
    for i = 1:cv.NumTestSets
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        
        % << 修改：使用适合多分类的模型，如 fitcecoc (推荐) 或 fitmnr >>
        % fitcecoc (Error-Correcting Output Codes) 是一个通用的多分类框架
        template = templateSVM('KernelFunction','linear','Standardize',true);
        mdl = fitcecoc(X(trainIdx,:), Y(trainIdx), ...
                       'Learners', template, ...
                       'Coding', 'onevsone');  % 多分类
        
        Y_pred = predict(mdl, X(testIdx,:));
        acc_fold = mean(Y_pred == Y(testIdx));
        accuracy_sum = accuracy_sum + acc_fold;
        confMat = confusionmat(Y(testIdx), Y_pred, ...
                               'Order', unique(Y)); 
        confMat_total = confMat_total + confMat;
    end
    
    avg_accuracy = accuracy_sum / cv.NumTestSets;
    fprintf('Features: %-20s | Average Accuracy: %.2f%%\n', feature_name{1}, avg_accuracy * 100);
        figure;
    cm = confusionchart(confMat_total, unique(Y));
    cm.Title = sprintf('Confusion Matrix - %s', feature_name{1});
    cm.RowSummary = 'row-normalized';
    cm.ColumnSummary = 'column-normalized';
end

fprintf('============================================================\n');
disp('To explore other models, type: classificationLearner(all_features_matrix, Y)');

% --- 辅助函数：计算互信息 ---
% 需要将此函数放在脚本末尾，或者保存为单独的 mi.m 文件
function I = mi(A, B)
    if size(A) ~= size(B)
        error('A and B must be of the same size.');
    end
    A = A(:); B = B(:);
    
    N = length(A);
    
    % 联合概率
    joint_counts = accumarray([A, B], 1);
    joint_prob = joint_counts / N;
    
    % 边缘概率
    prob_A = sum(joint_prob, 2);
    prob_B = sum(joint_prob, 1);
    
    % 互信息
    [ia, ja] = find(joint_prob > 0);
    H_A = -sum(prob_A(prob_A > 0) .* log2(prob_A(prob_A > 0)));
    H_B = -sum(prob_B(prob_B > 0) .* log2(prob_B(prob_B > 0)));
    
    joint_prob_vec = joint_prob(joint_prob > 0);
    H_AB = -sum(joint_prob_vec .* log2(joint_prob_vec));
    
    I = H_A + H_B - H_AB;
end