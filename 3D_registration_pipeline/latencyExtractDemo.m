%This program is used for plotting the lick latency change in a trial and across trials

%% 修改程序以支持选择多个文件并分别存储metaLickLatency变量
% 1. load the selected discrimination training data
[fileNames, path] = uigetfile('*.mat', 'Select the data files', 'MultiSelect', 'on');

if ischar(fileNames)
    fileNames = {fileNames}; % 如果只选择了一个文件，将其转换为单元数组
end

metaLickLatencyAll = cell(1, numel(fileNames)); % 创建一个单元数组存储每个文件的metaLickLatency

for fileIdx = 1:numel(fileNames)
    fileName = fileNames{fileIdx};
    fullPath = fullfile(path, fileName);
    
    % 加载数据文件
    data = load(fullPath);
    h = data.h; % 假设数据文件中包含变量h

    indexMatrix = h.data1(:,2);

    for i = 1:numel(indexMatrix)
        if mod(indexMatrix(i), 2) == 0
            indexMatrix(i) = 0;
        end
    end 

    indexMatrix = nonzeros(indexMatrix);

    % create the matrix to store the latency change only in Hit trial
    metalickLatency = zeros(length(indexMatrix), 2);
    metalickLatency(:, 1) = indexMatrix;
    for i = 1:length(indexMatrix)
        if indexMatrix(i) == 1 % only get the first lick latency in Hit trial
            metalickLatency(i, 2) = h.lickdata{i}(1, 2);
            if metalickLatency(i, 2) > 10
                metalickLatency(i, 2) = 0.005;
            end
        end
    end
    hitLickLatency = nonzeros(metalickLatency(:,2));
    metaLickLatencyAll{fileIdx} = hitLickLatency; % 存储当前文件的metaLickLatency
end

% 现在metaLickLatencyAll包含了每个文件的metaLickLatency

%% 修改散点图为带有钟形抖动的样式
figure;
hold on;
grayColor = [0.8, 0.8, 0.8]; % 灰色

% 添加条形图以指示每一天的中位数和标准差
for day = 1:10
    if day <= length(metaLickLatencyAll)
        lickLatency = metaLickLatencyAll{day}(all(metaLickLatencyAll{day}, 2), :); % 移除零值列
        x = day * ones(size(lickLatency)); % x轴为天数
        y = lickLatency; % y轴为lick latency
        jitter = 0.05 * randn(size(x)); % 钟形抖动
        scatter(x + jitter, y, 36, grayColor, 'filled'); % 绘制带抖动的灰色散点图

        % 计算中位数和标准差
        medianValue = median(y);
        stdValue = std(y);

        % 绘制中位数条形图
        line([day - 0.2, day + 0.2], [medianValue, medianValue], 'Color', 'k', 'LineWidth', 3);

        % 绘制标准差条形图
        line([day, day], [medianValue - stdValue, medianValue + stdValue], 'Color', 'k', 'LineWidth', 1);
    end
end

% 在图像顶部标注每一天的Hit试验数量
for day = 1:10
    if day <= length(metaLickLatencyAll)
        lickLatency = metaLickLatencyAll{day}(all(metaLickLatencyAll{day}, 2), :); % 移除零值列
        hitTrialCount = size(lickLatency, 1); % 计算Hit试验数量

        % 在图像顶部标注Hit试验数量
        text(day, 3.7, sprintf('%d', hitTrialCount), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', 'k');
    end
end

xlabel('Days');
ylabel('Lick Latency (s)');
title('Daily Lick Latency Scatter Plot with Bell-Shape Jitter');
set(gca, 'TickLength', [0.005, 0.005]);
xlim([0 11]);
ylim([-1 4])
xticks(0:10);
hold off;

%% 手动选择每只动物的.mat文件并以文件名标注在图上
figure;
hold on;
colors = lines(10); % 使用不同颜色表示四只动物

[fileNames, path] = uigetfile('*.mat', 'Select the data files for each animal', 'MultiSelect', 'on');
if ischar(fileNames)
    fileNames = {fileNames}; % 如果只选择了一个文件，将其转换为单元数组
end

for animalIdx = 1:length(fileNames)
    fileName = fileNames{animalIdx}; %
    fullPath = fullfile(path, fileName);
    data = load(fullPath);

    % 假设数据结构与之前一致
    [~, truncatedFileName, ~] = fileparts(fileName); 
    tempFieldname = fieldnames(data);
    hitLickLatencyAll.(truncatedFileName) = data.(tempFieldname{1});

    medianValues = zeros(1, 10);
    stdValues = zeros(1, 10);
    for day = 1:10
        if day <= length(hitLickLatencyAll.(truncatedFileName)) % 检查数据是否存在
            tempLatency = hitLickLatencyAll.(truncatedFileName){day};
            medianValues(day) = median(tempLatency); % 计算中位数
            stdValues(day) = std(tempLatency); % 计算标准差
        else
            medianValues(day) = NaN; % 如果数据不存在，用NaN填充
            stdValues(day) = NaN;
        end
    end
    
    if startsWith(fileName, 'WT')
        plot(1:10, medianValues, 'o-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 4, 'DisplayName', fileName); % Blue for WT
    elseif startsWith(fileName, 'APP')
        plot(1:10, medianValues, 'o-', 'Color', colors(animalIdx, :), 'LineWidth', 4, 'DisplayName', fileName); % Red for APP
    else
        plot(1:10, medianValues, 'o-', 'Color', colors(animalIdx, :), 'LineWidth', 4, 'DisplayName', fileName); % Default color
    end
    
    % 绘制标准差范围为误差棒样式
    %errorbar(1:10, medianValues, stdValues, 'Color', colors(animalIdx, :), 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
end

xlabel('Days');
ylabel('Median Lick Latency (s)');
title('Median Lick Latency Across Days for Selected Animals');
legend('show');
hold off;





