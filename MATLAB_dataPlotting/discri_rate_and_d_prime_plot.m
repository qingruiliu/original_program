


function discri_rate_and_d_prime_plot()
    % 使用GUI选择多个文件
    [fileNames, filePath] = uigetfile('*.mat', '选择实验数据文件', 'MultiSelect', 'on');
    
    % 处理用户取消选择的情况
    if isequal(fileNames, 0)
        disp('用户取消了文件选择');
        return;
    end
    
    % 确保fileNames是cell数组
    if ~iscell(fileNames)
        fileNames = {fileNames};
    end
    
    % 按文件名排序（假设文件名包含日期信息）
    fileNames = sort(fileNames);
    numDays = length(fileNames);
    
    % 初始化数据结构
    experimentData = struct('fileName', {}, 'date', {}, 'Nhit', [], ...
                           'Nmiss', [], 'Nfa', [], 'Ncr', [], ...
                           'correct_rate', [], 'd_prime', [], 'Phit', [], 'Pfa', []);
    
    % 处理每个文件
    for i = 1:numDays
        fullFilePath = fullfile(filePath, fileNames{i});
        fprintf('正在处理文件 %d/%d: %s\n', i, numDays, fileNames{i});
        load(fullFilePath);
        
        % 提取试验结果
        trialResult = h.data1(:,2);
        Nhit = sum(trialResult == 1);
        Nmiss = sum(trialResult == 2);
        Nfa = sum(trialResult == 3);
        Ncr = sum(trialResult == 4);
        correct_rate = (Nhit + Ncr) / length(trialResult);
        
        % 计算 d'
        Phit = Nhit / (Nhit + Nmiss);
        Pfa = Nfa / (Nfa + Ncr);
        
        % 处理边界情况（避免inf值）
        if Phit == 1
            Phit = 1 - 1/(2*(Nhit + Nmiss));
        elseif Phit == 0
            Phit = 1/(2*(Nhit + Nmiss));
        end
        
        if Pfa == 1
            Pfa = 1 - 1/(2*(Nfa + Ncr));
        elseif Pfa == 0
            Pfa = 1/(2*(Nfa + Ncr));
        end
        
        d_prime = norminv(Phit) - norminv(Pfa);
        
        % 保存到数据结构
        experimentData(i).fileName = fileNames{i};
        experimentData(i).date = i;  % 或从文件名提取日期
        experimentData(i).Nhit = Nhit;
        experimentData(i).Nmiss = Nmiss;
        experimentData(i).Nfa = Nfa;
        experimentData(i).Ncr = Ncr;
        experimentData(i).correct_rate = correct_rate;
        experimentData(i).d_prime = d_prime;
        experimentData(i).Phit = Phit;
        experimentData(i).Pfa = Pfa;
    end
    
    % 保存数据结构
    saveFileName = fullfile(filePath, 'multi_day_analysis.mat');
    save(saveFileName, 'experimentData');
    fprintf('\n已保存数据到: %s\n', saveFileName);
    
    % 绘制结果
    plotResults(experimentData, filePath);
    
    % 显示统计信息
    displayStatistics(experimentData);
end

function plotResults(experimentData, savePath)
    % 提取数据用于绘图
    days = [experimentData.date];
    correct_rates = [experimentData.correct_rate];
    d_primes = [experimentData.d_prime];
    
    % 绘制正确率图
    fig1 = figure('Name', 'Correct Rate Over Time', 'Position', [100, 100, 600, 500]);
    plot(days, correct_rates * 100, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.4470 0.7410]);
    hold on;
    plot(days, correct_rates * 100, 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 4);
    xlabel('Day (#)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Correct rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Correct Rate Over Time', 'FontSize', 14, 'FontWeight', 'bold');
    grid off;
    ylim([45, 100]);
    xlim([0.5, max(days)+0.5]);
    
    % 保存正确率图
    correctRateFileName = fullfile(savePath, 'correct_rate_over_time.png');
    saveas(fig1, correctRateFileName);
    fprintf('正确率图已保存为: %s\n', correctRateFileName);
    
    % 绘制 d' 图
    fig2 = figure('Name', 'd'' Over Time', 'Position', [750, 100, 600, 500]);
    plot(days, d_primes, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8500 0.3250 0.0980]);
    hold on;
    plot(days, d_primes, 'o', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 4);
    xlabel('Day (#)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('d''', 'FontSize', 12, 'FontWeight', 'bold');
    title('d'' Over Time', 'FontSize', 14, 'FontWeight', 'bold');
    grid off;
    xlim([0.5, max(days)+0.5]);
    ylim([-0.5 , max(d_primes)+0.5]);
    
    % 保存 d' 图
    dPrimeFileName = fullfile(savePath, 'd_prime_over_time.png');
    saveas(fig2, dPrimeFileName);
    fprintf('d''图已保存为: %s\n', dPrimeFileName);
end

function displayStatistics(experimentData)
    fprintf('\n========== 统计摘要 ==========\n');
    fprintf('总实验天数: %d\n', length(experimentData));
    fprintf('平均正确率: %.2f%%\n', mean([experimentData.correct_rate]) * 100);
    fprintf('正确率范围: %.2f%% - %.2f%%\n', ...
        min([experimentData.correct_rate]) * 100, max([experimentData.correct_rate]) * 100);
    fprintf('正确率标准差: %.2f%%\n', std([experimentData.correct_rate]) * 100);
    fprintf('\n');
    fprintf('平均 d'': %.2f\n', mean([experimentData.d_prime]));
    fprintf('d''范围: %.2f - %.2f\n', min([experimentData.d_prime]), max([experimentData.d_prime]));
    fprintf('d''标准差: %.2f\n', std([experimentData.d_prime]));
    fprintf('==============================\n\n');
    
    % 显示每天的详细数据
    fprintf('每天详细数据:\n');
    fprintf('%-20s\t正确率\t\td''\n', '文件名');
    fprintf('-----------------------------------------------------------\n');
    for i = 1:length(experimentData)
        fprintf('%-20s\t%.2f%%\t\t%.2f\n', ...
            experimentData(i).fileName, ...
            experimentData(i).correct_rate * 100, ...
            experimentData(i).d_prime);
    end
end