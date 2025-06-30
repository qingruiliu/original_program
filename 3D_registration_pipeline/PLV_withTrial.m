%% choose the demo .mat file
[fileName, path] = uigetfile('*.mat', 'Select the demo data for analysis', 'MultiSelect', 'off');
if isequal(fileName, 0)
    disp('User canceled the file selection.');
    return;
end
cd(path);
% Load the selected .mat file
load(fileName);

% the variable is "ROISegTraceTable"
cell1 = ROISegTraceTable.Segmented_trace{1};
cell2 = ROISegTraceTable.Segmented_trace{2};

% get the traces of the first 20 trials for each cell
cell1_trace = cell1.traces(1:20);
cell2_trace = cell2.traces(1:20);

% in each cell, 2nd column saves the time stamp of each trial, 3rd column saves the z-scored fluorescence

%create the cell array to store the interpolated traces
cell1_interp = cell(length(cell1_trace), 1);
cell2_interp = cell(length(cell2_trace), 1);

% interpolate the traces to have the same length, upsample to 10Hz
fs = 10;
dt = 1/fs;
% 可视化插值后的trace，并将每个trial的trace绘制在5x2子图中
num_trials = length(cell1_trace);
num_rows = 5;
num_cols = 2;

% 支持num_trials超过10时自动分页绘图
trials_per_fig = num_rows * num_cols;
num_figs = ceil(num_trials / trials_per_fig);
for fig_idx = 1:num_figs
    figure;
    for sub_idx = 1:trials_per_fig
        trial_idx = (fig_idx-1)*trials_per_fig + sub_idx;
        if trial_idx > num_trials
            break;
        end
        % 获取两个cell的最小和最大时间戳
        min_time = min([cell1_trace{trial_idx}(:, 2); cell2_trace{trial_idx}(:, 2)]);
        max_time = max([cell1_trace{trial_idx}(:, 2); cell2_trace{trial_idx}(:, 2)]);
        % 创建10Hz采样的公共时间轴
        common_axis = min_time:dt:max_time;
        % 对两个cell的trace进行线性插值
        cell1_interp{trial_idx} = [common_axis', interp1(cell1_trace{trial_idx}(:, 2), cell1_trace{trial_idx}(:, 3), common_axis, 'linear', 'extrap')'];
        cell2_interp{trial_idx} = [common_axis', interp1(cell2_trace{trial_idx}(:, 2), cell2_trace{trial_idx}(:, 3), common_axis, 'linear', 'extrap')'];
        % 绘制每个trial的trace到对应子图
        subplot(num_rows, num_cols, sub_idx);
        plot(cell1_interp{trial_idx}(:,1), cell1_interp{trial_idx}(:,2), 'b', 'LineWidth', 1.2); hold on;
        plot(cell2_interp{trial_idx}(:,1), cell2_interp{trial_idx}(:,2), 'r', 'LineWidth', 1.2);
        xlabel('Time (s)'); ylabel('Z-score');
        title(['Trial ' num2str(trial_idx)]);
        xregion([0 1]); % Set x-axis limits to -1s to 5s
        xlim([-1 8]);
        legend({'Cell1','Cell2'});
        grid off;
    end
    sgtitle(['Interpolated traces of each trial (Cell1: blue, Cell2: red), Fig ' num2str(fig_idx)]);
end
sgtitle('Interpolated traces of each trial (Cell1: blue, Cell2: red)');

%prepare the bandpassed filter parameters
low_cutoff = 0.1; % low cutoff frequency in Hz
high_cutoff = 1.15; % high cutoff frequency in Hz

% design a Butterworth bandpass filter
[b, a] = butter(2, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

% apply the filter 
for i = 1:length(cell1_interp)
    % Filter the traces
    cell1_interp{i}(:, 3) = filtfilt(b, a, cell1_interp{i}(:, 2));
    cell2_interp{i}(:, 3) = filtfilt(b, a, cell2_interp{i}(:, 2));

end
% calculate the instantaneous phase during the analysis window (-1s to 5s) of each trial
% using the Hilbert transform
cell1_phase = cell(length(cell1_interp), 1);
cell2_phase = cell(length(cell2_interp), 1);

% 在分析窗口内绘制滤波后的神经活动
figure;
for i = 1:length(cell1_interp)
    % 提取分析窗口 (0s 到 5s)
    analysis_window = cell1_interp{i}(:, 1) >= 0 & cell1_interp{i}(:, 1) <= 5;
    % 绘制滤波后的trace
    subplot(5, 4, i); % 每页最多20个trial
    plot(cell1_interp{i}(analysis_window, 1), cell1_interp{i}(analysis_window, 3), 'b', 'LineWidth', 1.2); hold on;
    plot(cell2_interp{i}(analysis_window, 1), cell2_interp{i}(analysis_window, 3), 'r', 'LineWidth', 1.2);
    xlabel('Time (s)'); ylabel('Filtered Z-score');
    title(['Filtered Trial ' num2str(i)]);
    legend({'Cell1','Cell2'});
    grid off;
    % Hilbert变换获取瞬时相位
    cell1_phase{i} = angle(hilbert(cell1_interp{i}(analysis_window, 3)));
    cell2_phase{i} = angle(hilbert(cell2_interp{i}(analysis_window, 3)));
end
sgtitle('Filtered neural activity in analysis window (0~5s)');
%calculate the Phase Locking Value (PLV) in each trial
PLV = zeros(length(cell1_interp), 1);
for i = 1:length(cell1_interp)
    % Calculate the phase difference between the two cells
    phase_diff = cell1_phase{i} - cell2_phase{i};
    % Calculate the PLV
    PLV(i) = abs(mean(exp(1i * phase_diff))); % Mean of the complex exponentials
end

%% plot the PLV change across trials, and label the trial contrast and result
figure;
plot(PLV, 'ok-','LineWidth', 1.5); % Plot PLV with black circles and lines
xlabel('Trial Number');
ylabel('Phase Locking Value (PLV)');
%title('PLV between Cell 1 and Cell 2 across Trials');
xlim([0 21]);

% Add trial contrast and result labels
result = cell1.trialResult(1:20); % Assuming trialResult is a field in the cell structure
contrast = cell1.trialContrast(1:20); % Assuming trialContrast is a field in the cell structure
%plot the trial result as the colored filled area, and label the contrast by text
hold on;
hitColor = [0.25 0.8 0.25]; % Color for hit trials
missColor = [1 0.54 0.1]; % Color for miss trials
falseAlarmColor = [0.83 0.14 0.14]; % Color for false alarm trials
correctRejectColor = [0.27 0.25 0.8]; % Color for correct rejection trials
for i = 1:length(result)
    if result(i) == 1 % Hit trial
        fill([i-0.05, i+0.05, i+0.05, i-0.05], [0, 0, max(PLV), max(PLV)], hitColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    elseif result(i) == 2 % Miss trial
        fill([i-0.05, i+0.05, i+0.05, i-0.05], [0, 0, max(PLV), max(PLV)], missColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    elseif result(i) == 3 % False alarm trial
        fill([i-0.05, i+0.05, i+0.05, i-0.05], [0, 0, max(PLV), max(PLV)], falseAlarmColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    elseif result(i) == 4 % Correct rejection trial
        fill([i-0.05, i+0.05, i+0.05, i-0.05], [0, 0, max(PLV), max(PLV)], correctRejectColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    % Add text label for contrast
    text(i, max(PLV) + 0.05, num2str(contrast(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
hold off;   
