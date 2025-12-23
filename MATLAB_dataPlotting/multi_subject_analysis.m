function multi_subject_analysis()
    % 读取多个实验个体的数据并汇总在同一张表格上
    
    % 使用GUI选择多个个体的数据文件夹或文件
    fprintf('请选择第一个个体的multi_day_analysis.mat文件\n');
    
    % 初始化数据存储
    allSubjectsData = {};
    subjectNames = {};
    continueSelection = true;
    subjectCount = 0;
    
    % 循环选择多个个体的数据
    while continueSelection
        [fileName, filePath] = uigetfile('*.mat', ...
            sprintf('选择个体 %d 的 multi_day_analysis.mat 文件 (取消结束选择)', subjectCount + 1));
        
        if isequal(fileName, 0)
            % 用户取消选择
            if subjectCount == 0
                disp('未选择任何文件，程序结束');
                return;
            else
                break; % 结束选择
            end
        end
        
        subjectCount = subjectCount + 1;
        
        % 加载数据
        fullFilePath = fullfile(filePath, fileName);
        data = load(fullFilePath);
        
        % 提取个体名称（从文件路径或让用户输入）
        prompt = sprintf('请输入个体 %d 的名称 (默认: Subject_%d): ', subjectCount, subjectCount);
        subjectName = input(prompt, 's');
        if isempty(subjectName)
            subjectName = sprintf('Subject_%d', subjectCount);
        end
        
        % 保存数据
        allSubjectsData{subjectCount} = data.experimentData;
        subjectNames{subjectCount} = subjectName;
        
        fprintf('已加载个体: %s (%d天数据)\n', subjectName, length(data.experimentData));
    end
    
    % 创建汇总表格
    summaryTable = createSummaryTable(allSubjectsData, subjectNames);
    
    % 显示表格
    disp('========== 多个体数据汇总 ==========');
    disp(summaryTable);
    
    % 保存表格
    [saveFileName, savePath] = uiputfile('*.xlsx', '保存汇总表格', 'multi_subject_summary.xlsx');
    if ~isequal(saveFileName, 0)
        fullSavePath = fullfile(savePath, saveFileName);
        writetable(summaryTable, fullSavePath);
        fprintf('\n表格已保存为: %s\n', fullSavePath);
    end
    
    % 绘制对比图
    plotComparison(allSubjectsData, subjectNames, savePath);
end

function summaryTable = createSummaryTable(allSubjectsData, subjectNames)
    % 创建汇总表格
    
    numSubjects = length(allSubjectsData);
    
    % 初始化表格数据
    subjectNameCol = {};
    numDaysCol = [];
    avgCorrectRateCol = [];
    stdCorrectRateCol = [];
    minCorrectRateCol = [];
    maxCorrectRateCol = [];
    avgDPrimeCol = [];
    stdDPrimeCol = [];
    minDPrimeCol = [];
    maxDPrimeCol = [];
    
    % 计算每个个体的统计数据
    for i = 1:numSubjects
        data = allSubjectsData{i};
        
        correct_rates = [data.correct_rate];
        d_primes = [data.d_prime];
        
        subjectNameCol{i} = subjectNames{i};
        numDaysCol(i) = length(data);
        avgCorrectRateCol(i) = mean(correct_rates) * 100;
        stdCorrectRateCol(i) = std(correct_rates) * 100;
        minCorrectRateCol(i) = min(correct_rates) * 100;
        maxCorrectRateCol(i) = max(correct_rates) * 100;
        avgDPrimeCol(i) = mean(d_primes);
        stdDPrimeCol(i) = std(d_primes);
        minDPrimeCol(i) = min(d_primes);
        maxDPrimeCol(i) = max(d_primes);
    end
    
    % 创建表格
    summaryTable = table(subjectNameCol', numDaysCol', ...
        avgCorrectRateCol', stdCorrectRateCol', minCorrectRateCol', maxCorrectRateCol', ...
        avgDPrimeCol', stdDPrimeCol', minDPrimeCol', maxDPrimeCol', ...
        'VariableNames', {'Subject', 'NumDays', ...
        'AvgCorrectRate', 'StdCorrectRate', 'MinCorrectRate', 'MaxCorrectRate', ...
        'AvgDPrime', 'StdDPrime', 'MinDPrime', 'MaxDPrime'});
end

function plotComparison(allSubjectsData, subjectNames, savePath)
    % 绘制多个体对比图
    
    numSubjects = length(allSubjectsData);
    colors = lines(numSubjects); % 自动生成不同颜色
    
    % 绘制正确率对比图
    fig1 = figure('Name', 'Multi-Subject Correct Rate Comparison', 'Position', [100, 100, 800, 600]);
    hold on;
    
    for i = 1:numSubjects
        data = allSubjectsData{i};
        days = [data.date];
        correct_rates = [data.correct_rate] * 100;
        
        plot(days, correct_rates, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
            'Color', colors(i,:), 'DisplayName', subjectNames{i});
    end
    
    xlabel('Day (#)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Correct Rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Multi-Subject Correct Rate Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    %legend('Location', 'best');
    grid on;
    ylim([40, 100]);
    hold off;
    
    % 保存正确率对比图
    if ~isequal(savePath, 0)
        correctRateFileName = fullfile(savePath, 'multi_subject_correct_rate.png');
        saveas(fig1, correctRateFileName);
        fprintf('多个体正确率对比图已保存为: %s\n', correctRateFileName);
    end
    
    % 绘制 d' 对比图
    fig2 = figure('Name', 'Multi-Subject d'' Comparison', 'Position', [150, 150, 800, 600]);
    hold on;
    
    for i = 1:numSubjects
        data = allSubjectsData{i};
        days = [data.date];
        d_primes = [data.d_prime];
        
        plot(days, d_primes, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
            'Color', colors(i,:), 'DisplayName', subjectNames{i});
    end
    
    xlabel('Day (#)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('d''', 'FontSize', 12, 'FontWeight', 'bold');
    title('Multi-Subject d'' Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    %legend('Location', 'best');
    grid off;
    hold off;
    
    % 保存 d' 对比图
    if ~isequal(savePath, 0)
        dPrimeFileName = fullfile(savePath, 'multi_subject_d_prime.png');
        saveas(fig2, dPrimeFileName);
        fprintf('多个体d''对比图已保存为: %s\n', dPrimeFileName);
    end
    
    % 绘制箱线图对比
    fig3 = figure('Name', 'Multi-Subject Box Plot Comparison', 'Position', [200, 200, 1000, 500]);
    
    % 正确率箱线图
    subplot(1, 2, 1);
    correctRateData = [];
    groupLabels = [];
    for i = 1:numSubjects
        data = allSubjectsData{i};
        correct_rates = [data.correct_rate] * 100;
        correctRateData = [correctRateData, correct_rates];
        groupLabels = [groupLabels, repmat({subjectNames{i}}, 1, length(correct_rates))];
    end
    boxplot(correctRateData, groupLabels);
    ylabel('Correct Rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Correct Rate Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid off;
    
    % d' 箱线图
    subplot(1, 2, 2);
    dPrimeData = [];
    groupLabels = [];
    for i = 1:numSubjects
        data = allSubjectsData{i};
        d_primes = [data.d_prime];
        dPrimeData = [dPrimeData, d_primes];
        groupLabels = [groupLabels, repmat({subjectNames{i}}, 1, length(d_primes))];
    end
    boxplot(dPrimeData, groupLabels);
    ylabel('d''', 'FontSize', 12, 'FontWeight', 'bold');
    title('d'' Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid off;
    
    % 保存箱线图
    if ~isequal(savePath, 0)
        boxPlotFileName = fullfile(savePath, 'multi_subject_boxplot.png');
        saveas(fig3, boxPlotFileName);
        fprintf('多个体箱线图已保存为: %s\n', boxPlotFileName);
    end
end
