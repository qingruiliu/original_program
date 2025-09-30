%% choose the demo .mat file
[fileName, path] = uigetfile('*.mat', 'Select the demo data for analysis', 'MultiSelect', 'off');
if isequal(fileName, 0)
    disp('User canceled the file selection.');
    return;
end
cd(path);
% Load the selected .mat file
load(fileName);

% save the actual center and the transformed center of each ROI
center = ROISegTraceTable.Center;  %center is in Nx3 format (x,y,z), where N is the number of ROIs
transCenter = ROISegTraceTable.TransformedCenter; %transCenter is in Nx3 format (x,y,z), where N is the number of ROIs
idx = ROISegTraceTable.ROIIndex; % ROIIndex is a vector of indices corresponding to each ROI


% 计算每对ROI的距离，并保存配对索引、欧氏距离和切向距离到表格
pair_idx1 = [];
pair_idx2 = [];
pair_dist = [];
pair_tanDist = [];
for i = 1:size(transCenter,1)
    for j = i+1:size(transCenter,1)
        dist = norm(center(i,:) - center(j,:));
        tanDist = norm(transCenter(i,1:2) - transCenter(j,1:2));
        if tanDist <= 10 && dist > 20   % 只保留切向距离小于20且欧氏距离大于20的ROI对
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
%add another column for saving the pearson R
aligned_pairs_table.R = cell(height(aligned_pairs_table), 1);
%add another column for saving the PLV results in cell array
aligned_pairs_table.Features = cell(height(aligned_pairs_table), 1);



% 针对每一对aligned_pairs_table中的ROI对，进行分析并保存结果
num_pairs = height(aligned_pairs_table);
pair_PLV_table = table();
fs = 10;
dt = 1/fs;
low_cutoff = 0.1;
high_cutoff = 1.15;
[b, a] = butter(2, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

% 定义绘图函数
function plot_trials(trial_indices, condition_name, contrast_val_str, roi1, roi2, cell1_interp, cell2_interp)
    if isempty(trial_indices)
        fprintf('No trials found for condition: %s\n', condition_name);
        return;
    end
    num_trials_to_plot = length(trial_indices);
    figure('Name', sprintf('ROI%d-ROI%d Traces: %s', roi1, roi2, condition_name));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    for i = 1:num_trials_to_plot
        k = trial_indices(i);
        subplot(ceil(sqrt(num_trials_to_plot)), ceil(sqrt(num_trials_to_plot)), i);
        set(gca, 'FontSize', 12);
        plot(cell1_interp{k}(:,1), cell1_interp{k}(:,2), 'b', 'LineWidth', 1.2); hold on;
        plot(cell2_interp{k}(:,1), cell2_interp{k}(:,2), 'r', 'LineWidth', 1.2);
        xlabel('Time (s)', 'FontSize', 14); ylabel('Z-score', 'FontSize', 14);
        title(['Trial ' num2str(k)], 'FontSize', 14);
        xlim([-1 8]);
        xregion([0 1]); %time window of visual stimulation
        grid off;
    end
    sgtitle(sprintf('ROI%d-ROI%d: %s (Contrast: %s)', roi1, roi2, condition_name, contrast_val_str), 'FontSize', 20);
    uiwait(gcf);
end

wb = waitbar(0, 'Processing PLV for each pair of ROIs...');
% 遍历每一对ROI
for p = 1:num_pairs
    waitbar(p/num_pairs, wb, sprintf('Processing PLV for pair %d of %d', p, num_pairs));
    roi1 = aligned_pairs_table.ROI1(p);
    roi2 = aligned_pairs_table.ROI2(p);
    % 获取两个ROI的trace结构
    cell1 = ROISegTraceTable.Segmented_trace{find(ROISegTraceTable.ROIIndex == roi1,1)};
    cell2 = ROISegTraceTable.Segmented_trace{find(ROISegTraceTable.ROIIndex == roi2,1)};
    % 加入确定signal contamination的功能
    cell1_totalTrace = vertcat(cell1.traces{:}); 
    cell1_totalTrace = cell1_totalTrace(:, 3); % 取Z-score
    cell2_totalTrace = vertcat(cell2.traces{:});
    cell2_totalTrace = cell2_totalTrace(:, 3); % 取Z-score

    % 保持两个totalTrace为相同长度
    min_length = min(length(cell1_totalTrace), length(cell2_totalTrace));
    cell1_totalTrace = cell1_totalTrace(1:min_length);
    cell2_totalTrace = cell2_totalTrace(1:min_length);
    % 计算两个ROI的相关性
    tempR = corrcoef(cell1_totalTrace, cell2_totalTrace);
    tempR = tempR(1,2);

    %plot the scatter plot of two traces
    figure('Name','Signal Contamination Check');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'FontSize', 20);
    scatter(cell1_totalTrace, cell2_totalTrace, 36, 'filled', 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor', 'none');
    xlabel(sprintf('ROI%d Z-score', roi1), 'FontSize', 20);
    ylabel(sprintf('ROI%d Z-score', roi2), 'FontSize', 20);
    title(sprintf('Scatter Plot of ROI%d and ROI%d, R = %s', roi1, roi2,num2str(tempR)), 'FontSize', 20);

    % manually assign the result of signal contamination
    out = questdlg('Is there signal contamination in this pair?', ...
        'Signal Contamination Check', 'Yes', 'No', 'No');
    if strcmp(out, 'Yes')
        aligned_pairs_table.R{p} = NaN; % 如果有信号污染，相关系数设为NaN
        continue; % 跳过后续处理
    elseif strcmp(out, 'No')
        aligned_pairs_table.R{p} = tempR; % 保存相关系数到表
    end
    uiwait(gcf);
    trialNum = height(cell1.traces); % 获取trial数量
    % 取前20个trial
    cell1_trace = cell1.traces(1:trialNum);
    cell2_trace = cell2.traces(1:trialNum);
    cell1_interp = cell(length(cell1_trace), 1);
    cell2_interp = cell(length(cell2_trace), 1);
    PLV = zeros(length(cell1_trace), 1);
    TrialCorr = zeros(length(cell1_trace), 1);   

    for i = 1:length(cell1_trace)
        min_time = min([cell1_trace{i}(:, 2); cell2_trace{i}(:, 2)]);
        max_time = max([cell1_trace{i}(:, 2); cell2_trace{i}(:, 2)]);
        common_axis = min_time:dt:max_time;
        cell1_interp{i} = [common_axis', interp1(cell1_trace{i}(:, 2), cell1_trace{i}(:, 3), common_axis, 'linear', 'extrap')'];
        cell2_interp{i} = [common_axis', interp1(cell2_trace{i}(:, 2), cell2_trace{i}(:, 3), common_axis, 'linear', 'extrap')'];

        %=== new variable ===%
        % 分析窗口
        analysis_window = cell1_interp{i}(:, 1) >= -0.5 & cell1_interp{i}(:, 1) <= 5;
        trace1_window = cell1_interp{i}(analysis_window, 2); % 使用未滤波的Z-score
        trace2_window = cell2_interp{i}(analysis_window, 2);
        if length(trace1_window) > 1 && length(trace2_window) > 1 && ~any(isnan(trace1_window)) && ~any(isnan(trace2_window))
            corr_matrix = corrcoef(trace1_window, trace2_window);
            TrialCorr(i) = corr_matrix(1, 2);
        else
            TrialCorr(i) = NaN;
        end

        % === old PLV method ===
        cell1_interp{i}(:, 3) = filtfilt(b, a, cell1_interp{i}(:, 2));
        cell2_interp{i}(:, 3) = filtfilt(b, a, cell2_interp{i}(:, 2));
        % Hilbert变换获取瞬时相位
        cell1_phase = angle(hilbert(cell1_interp{i}(analysis_window, 3)));
        cell2_phase = angle(hilbert(cell2_interp{i}(analysis_window, 3)));
        % 计算PLV
        phase_diff = cell1_phase - cell2_phase;
        PLV(i) = abs(mean(exp(1i * phase_diff)));
    end
    
    % 获取不同条件下的trial索引
    result = cell1.trialResult(1:length(PLV));
    contrast = cell1.trialContrast(1:length(PLV));
    
    idx_100_hit = find(result == 1 & contrast == 1);
    idx_10_hit = find(result == 1 & contrast == 0.1);
    idx_100_cr = find(result == 4 & contrast == 1); % CR at 100% contrast (FA)
    idx_10_cr = find(result == 4 & contrast == 0.1); % CR at 10% contrast (FA)

    %% 绘制四张图
    plot_trials(idx_100_hit, '100% Hit', '1', roi1, roi2, cell1_interp, cell2_interp);
    plot_trials(idx_10_hit, '10% Hit', '0.1', roi1, roi2, cell1_interp, cell2_interp);
    plot_trials(idx_100_cr, '100% CR', '1', roi1, roi2, cell1_interp, cell2_interp);
    plot_trials(idx_10_cr, '10% CR', '0.1', roi1, roi2, cell1_interp, cell2_interp);


    % 可视化PLV动态变化及trial结果
    figure('Name',sprintf('ROI%d-ROI%d PLV',roi1,roi2));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'FontSize', 20);
    plot(PLV, 'ok-','LineWidth', 1.5);
    xlabel('Trial Number', 'FontSize', 20);
    ylabel('Phase Locking Value (PLV)', 'FontSize', 20);
    title(sprintf('PLV between ROI%d and ROI%d across Trials',roi1,roi2), 'FontSize', 20);
    xlim([0 length(PLV)+1]);
    ylim([0 1]);
    box off
    % 叠加trial结果色条和对比度标签
    result = cell1.trialResult(1:length(PLV)); % 兼容trial数变化
    contrast = cell1.trialContrast(1:length(PLV));
    hold on;
    hitColor = [0.25 0.8 0.25];
    missColor = [1 0.54 0.1];
    falseAlarmColor = [0.83 0.14 0.14];
    correctRejectColor = [0.27 0.25 0.8];
    for iBar = 1:length(result)
        if result(iBar) == 1
            fill([iBar-0.05, iBar+0.05, iBar+0.05, iBar-0.05], [0, 0, max(PLV), max(PLV)], hitColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        elseif result(iBar) == 2
            fill([iBar-0.05, iBar+0.05, iBar+0.05, iBar-0.05], [0, 0, max(PLV), max(PLV)], missColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        elseif result(iBar) == 3
            fill([iBar-0.05, iBar+0.05, iBar+0.05, iBar-0.05], [0, 0, max(PLV), max(PLV)], falseAlarmColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        elseif result(iBar) == 4
            fill([iBar-0.05, iBar+0.05, iBar+0.05, iBar-0.05], [0, 0, max(PLV), max(PLV)], correctRejectColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        % Add text label for contrast
        text(iBar, max(PLV) + 0.05, num2str(contrast(iBar)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20);
    end

    %将160个trial的PLV结果保存在aligned_pairs_table中，形式为1x3，包括result，contrast和PLV
    PLV_data = [result,contrast,PLV];
    aligned_pairs_table.PLV{p} = PLV_data; % 保存PLV数据到表格中
    uiwait(gcf); % 等待用户关闭图形窗口
    
    % 绘制PLV随contrast变化的分组散点+均值线图（横轴为对数contrast，颜色区分trial结果）
    figure('Name',sprintf('ROI%d-ROI%d PLV by Contrast & Result',roi1,roi2));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'FontSize', 20);
    hold on;
    result_types = [1 2 3 4];
    result_labels = {'Hit','Miss','FA','CR'};
    contrast_values = unique(contrast);
    % 使用与前面一致的配色方案
    colors = [0.25 0.8 0.25; 1 0.54 0.1; 0.83 0.14 0.14; 0.27 0.25 0.8];
    markerstyles = {'o','s','^','d'};
    width = 0.15; % 横向微移宽度
    for r = 1:length(result_types)
        y_means = nan(size(contrast_values));
        for c = 1:length(contrast_values)
            idx = (result == result_types(r)) & (contrast == contrast_values(c));
            % 横向微移
            x = log10(contrast_values(c)) + (r-2.5)*width/2;
            scatter(repmat(x, sum(idx), 1), PLV(idx), 36, ...
                'MarkerFaceColor', colors(r,:), ...
                'MarkerEdgeColor', colors(r,:), ...
                'MarkerFaceAlpha', 0.7, ...
                'Marker', markerstyles{r});
            % 计算均值
            if sum(idx)>0
                y_means(c) = mean(PLV(idx));
            end
        end
        % 画均值线
        plot(log10(contrast_values) + (r-2.5)*width/2, y_means, '-', ...
            'Color', colors(r,:), ...
            'LineWidth',2, ...
            'Marker', markerstyles{r}, ...
            'MarkerFaceColor',colors(r,:), ...
            'MarkerSize', 8);
    end
    % 设置X轴为对数分布
    set(gca, 'XTick', log10(contrast_values), 'XTickLabel', arrayfun(@num2str, contrast_values, 'UniformOutput', false), 'FontSize', 20);
    xlabel('Contrast (log scale)', 'FontSize', 20);
    ylabel('PLV', 'FontSize', 20);
    title(sprintf('PLV vs Contrast (ROI%d-ROI%d, color=result)',roi1,roi2), 'FontSize', 20);
    % 调整图例的颜色
    h = zeros(1,4);
    for r = 1:4
        h(r) = plot(nan, nan, markerstyles{r}, ...
            'MarkerFaceColor', colors(r,:), ...
            'MarkerEdgeColor', colors(r,:), ...
            'Color', colors(r,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 8);
    end
    legend(h, result_labels, 'Location', 'best', 'FontSize', 20);
    box off; grid on;
    hold off;
    uiwait(gcf); % 等待用户关闭图形窗口

    % << 修改：将所有特征保存到一个结构化的表中 >>
    feature_data = table(result, contrast, PLV, TrialCorr, 'VariableNames', {'Result', 'Contrast', 'PLV', 'TrialCorr'});
    aligned_pairs_table.Features{p} = feature_data; % 保存特征数据到表格中

end
close(wb);
disp('Feature extraction completed.');

%% Part 2: Decoding Analysis
disp('Starting decoding analysis...');

% 1. 准备数据集
% 从 aligned_pairs_table 中提取所有 trial 的特征和标签
all_features = [];
all_labels_decision = [];
all_labels_contrast = [];
all_labels_correctness = [];

for p = 1:height(aligned_pairs_table)
    % 跳过有信号污染或没有数据的对
    if isempty(aligned_pairs_table.Features{p})
        continue;
    end
    
    pair_features_table = aligned_pairs_table.Features{p};
    
    % 移除包含NaN的行 (例如，相关性无法计算的trial)
    pair_features_table = rmmissing(pair_features_table);
    
    if isempty(pair_features_table)
        continue;
    end

    % 选择用于解码的特征。这里我们使用PLV和TrialCorr
    features_for_this_pair = [pair_features_table.PLV, pair_features_table.TrialCorr];
    
    % 获取标签
    decision_labels = pair_features_table.Result;
    contrast_labels = pair_features_table.Contrast;
    
    % 创建正确性标签 (1 for Correct, 0 for Incorrect)
    % Hit(1) 和 Correct Reject(4) 是正确的
    correctness_labels = (decision_labels == 1 | decision_labels == 4);

    % 汇总所有神经元对的数据
    all_features = [all_features; features_for_this_pair];
    all_labels_decision = [all_labels_decision; decision_labels];
    all_labels_contrast = [all_labels_contrast; contrast_labels];
    all_labels_correctness = [all_labels_correctness; correctness_labels];
end

if isempty(all_features)
    disp('No valid data available for decoding. Exiting.');
    return;
end

% 2. 训练和评估分类器
% 示例：解码动物的选择是否正确 (Correct vs. Incorrect)
disp('Training a model to decode choice correctness...');

% 特征矩阵 X 和 标签向量 Y
X = all_features; 
Y = all_labels_correctness;

% 检查是否有足够的数据和类别
if length(unique(Y)) < 2
    disp('Not enough classes in the data to perform classification.');
    return;
end

% 使用交叉验证来评估模型性能，防止过拟合
try
    cv = cvpartition(Y, 'KFold', 10); % 10折交叉验证
catch ME
    disp('Could not create cross-validation partition. Not enough data points?');
    disp(ME.message);
    return;
end

% 使用逻辑回归模型 (fitglm) 作为示例
accuracy_sum = 0;
for i = 1:cv.NumTestSets
    trainIdx = cv.training(i);
    testIdx = cv.test(i);
    
    % 训练模型
    mdl = fitglm(X(trainIdx,:), Y(trainIdx), 'Distribution', 'binomial', 'Link', 'logit');
    
    % 预测
    Y_pred_prob = predict(mdl, X(testIdx,:));
    Y_pred = Y_pred_prob > 0.5; % 将概率转换为类别 (0或1)
    
    % 计算准确率
    accuracy_sum = accuracy_sum + sum(Y_pred == Y(testIdx)) / length(Y(testIdx));
end

avg_accuracy = accuracy_sum / cv.NumTestSets;
fprintf('============================================================\n');
fprintf('Decoding Result (Correct vs. Incorrect Choice)\n');
fprintf('Features used: PLV and Trial-by-Trial Correlation\n');
fprintf('Model: Logistic Regression\n');
fprintf('Average decoding accuracy (10-fold CV): %.2f%%\n', avg_accuracy * 100);
fprintf('============================================================\n');

% 提示：您也可以使用MATLAB的 Classification Learner App 来快速尝试多种模型
% classificationLearner(X, Y);
disp('To explore other models, type: classificationLearner(X, Y)');
