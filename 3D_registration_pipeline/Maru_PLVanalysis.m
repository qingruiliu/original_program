%% pre-assign the data for analysis purpose
% choose the target data for analysis, multiple selection
[name,path] = uigetfile('*.mat','Select the data for analysis','MultiSelect','on');
cd(path)
fileNames = {'ROISegTraceTable_M3.mat', 'ROISegTraceTable_M9.mat', 'ROISegTraceTable_M17.mat'};

% 预分配结构体数组
numFiles = numel(fileNames);
dataStruct(numFiles) = struct('Data', [], 'OriginalFileName', '', 'SucessRate', []);

% 遍历每个文件并存入结构体
for i = 1:numFiles
    % 加载 .mat 文件中的数据
    fileData = load(fileNames{i});
    
    % 假设 .mat 文件中包含名为 'tableData' 的表格数据
    tableField = fieldnames(fileData);
    dataStruct(i).Data = fileData.(tableField{1}); % 这里假设第一个字段是所需数据
    
    % 存储文件名
    dataStruct(i).OriginalFileName = fileNames{i};

    behData = dataStruct(i).Data.Segmented_trace{1,1};
    behDataSort = sortrows(behData,{'trialContrast','trialResult'},{'ascend','ascend'});
    correctRate = zeros(1,4);
    for j = 1:4
        correctNum = behDataSort(behDataSort.trialContrast == 10^(j-4),:); 
        trialMark = correctNum.trialResult;
        correctRate(j) = sum(mod(trialMark,3) == 1);
    end
    dataStruct(i).SucessRate = correctRate/40;
end

% 查看结果
SumData = dataStruct;
%clearvars -except SumData


% pre-defined parameters
ClmDis=10;     % in-column pair distance [um]
AdClmDis=30;   % adjacent-column pair distance [um]
ExDis=20;      % threshold Euclidean distance [um]
Th=1;          % threshold for PLV
StartBin=1;    
EndBin=20;
Contrast=[0.001,0.01,0.1,1];
NTrial=size(SumData(1).Data.Segmented_trace{1,1},1); % number of trials
NAnimal=3;
NSurr=1000;

NTimePoints=EndBin;

% 時間軸
Fs = 2.37; % サンプリング周波数 [Hz]
time = (0:NTimePoints-1) / Fs; % time window (about 0~8s)


for iAnimal=1:3
    
    tic
    %細胞座標を処理しやすい形に変換
    SumData(iAnimal).CellPosiUM=cell2mat({SumData(iAnimal).Data.TransformedCenter});
    NCell=size(SumData(iAnimal).CellPosiUM,1);

    %Tangential方向の細胞間距離を算出
    TanDisTable=[];
    for iCell=1:NCell
        for jCell=1:NCell
            SumData(iAnimal).TanDisUMTable(iCell,jCell)=norm(SumData(iAnimal).CellPosiUM(iCell,1:2)-SumData(iAnimal).CellPosiUM(jCell,1:2));
        end
    end

    %直線距離を求める。
    NCell=size(SumData(iAnimal).CellPosiUM,1);
    for iCell=1:NCell
        for jCell=1:NCell
            SumData(iAnimal).StDisUMTable(iCell,jCell)=norm(SumData(iAnimal).CellPosiUM(iCell,:)-SumData(iAnimal).CellPosiUM(jCell,:));
        end
    end
    
    CaData_Trial=[];
    for iTrial=1:NTrial
        for iCell=1:NCell
            CaData_Trial(iCell,:,iTrial)=SumData(iAnimal).Data.Segmented_trace{iCell,1}.traces{iTrial,1}(StartBin:EndBin,3);
        end
    end


    PLV_Trial=zeros(NCell, NCell,NTrial);
    for iTrial=1:NTrial

        data=CaData_Trial(:,:,iTrial);

        % ----- フィルタリング -----
        % 周波数帯域を選択
        low_cutoff = 0.01;  % 下限 [Hz]
        high_cutoff = 1.15; % 上限 [Hz]
        [b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2)); % 2次バターワースフィルタ

        filtered_data = zeros(size(data));
        for i = 1:NCell
            filtered_data(i, :) = filtfilt(b, a, double(data(i, :))); % フィルタリング
        end

        % ----- ヒルベルト変換 -----
        % 瞬時位相の抽出
        instant_phase=[];
        for i = 1:NCell
            instant_phase(i,:) = angle(hilbert(filtered_data(i,:),20)); % 各細胞の瞬時位相
        end

        % ----- 位相ロッキング値 (PLV) の計算 -----
        % ペアごとのPLVを計算
        PLV = zeros(NCell, NCell); % PLV行列の初期化
        for i = 1:NCell
            for j = i+1:NCell
                % 位相差を計算
                phase_diff = instant_phase(i, :) - instant_phase(j, :);

                % 位相ロッキング値を計算
                PLV(i, j) = abs(mean(exp(1i * phase_diff))); % 平均して振幅を計算
                PLV(j, i) = PLV(i, j); % 対称性を持たせる
            end
        end

        PLV_Trial(:,:,iTrial)=PLV;
        SumData(iAnimal).PLV_Trial=PLV_Trial;

    end

    % ----- 可視化 -----
    % PLV行列をヒートマップで表示
    figure;
    imagesc(PLV);
    colorbar;
    title('Phase Locking Value (PLV)');
    xlabel('Neuron Index');
    ylabel('Neuron Index');
    axis square;

    m=1;
    n=1;
    ClmPLV_Trial=[];
    OtherPLV_Trial=[];
    TanDisTable=SumData(iAnimal).TanDisUMTable;
    StDisTable=SumData(iAnimal).StDisUMTable;
    for iCell=1:NCell
        for jCell=1:NCell
            if iCell~=jCell&iCell<jCell
                if TanDisTable(iCell,jCell)<=ClmDis&StDisTable(iCell,jCell)>=ExDis
                    ClmPLV_Trial(m,:)=squeeze(PLV_Trial(iCell,jCell,:));
                    m=m+1;
                end
                if TanDisTable(iCell,jCell)>AdClmDis&TanDisTable(iCell,jCell)<=(AdClmDis+ClmDis)&StDisTable(iCell,jCell)>=ExDis
                    OtherPLV_Trial(n,:)=squeeze(PLV_Trial(iCell,jCell,:));
                    n=n+1;
                end
            end
        end
        tocf(2)
    end
    
    SumData(iAnimal).ClmPLV_Trial=ClmPLV_Trial;
    SumData(iAnimal).OtherPLV_Trial=OtherPLV_Trial;
    
    %ResultとContrastのsequenceを抽出する。
    SumData(iAnimal).ResultSeq=SumData(iAnimal).Data.Segmented_trace{1,1}.trialResult;
    SumData(iAnimal).ContrastSeq=SumData(iAnimal).Data.Segmented_trace{1,1}.trialContrast;
    
    for iSurr=1:NSurr+1
        

        d2=SumData(iAnimal).ResultSeq;
        d3=SumData(iAnimal).ContrastSeq;

        for iResult=1:4
            for iCont=1:4
                TempIdx=[];
                TempIdx=find(d2==iResult&d3==Contrast(iCont));

                if isempty(TempIdx)==0
                    if iSurr==1
                        d=ClmPLV_Trial(:,TempIdx);
                    else
                        d=OtherPLV_Trial(randperm(n-1,m-1),TempIdx);
                    end
                    
                    TempCorr=[];
                    k=1;
                    for iTrial=1:length(TempIdx)
                        for jTrial=1:length(TempIdx)
                            if iTrial<jTrial
                                TempCorr(k)=corr(d(:,iTrial),d(:,jTrial));
                                k=k+1;
                            end
                        end
                    end
                    AveCorrBetTrial(iResult,iCont,iSurr)=mean(TempCorr);
                else
                    AveCorrBetTrial(iResult,iCont,iSurr)=NaN;
                end
            end
        end

    end
    
%     for iSurr=1:NSurr+1
%         
%         if iSurr==1
%             d2=SumData(iAnimal).ResultSeq;
%             d3=SumData(iAnimal).ContrastSeq;
%         else
%             d2_0=SumData(iAnimal).ResultSeq;
%             d2=d2_0(randperm(length(d2_0)));
%             d3_0=SumData(iAnimal).ContrastSeq;
%             d3=d3_0(randperm(length(d3_0)));
%         end
% 
%         for iResult=1:4
%             for iCont=1:4
%                 TempIdx=[];
%                 TempIdx=find(d2==iResult&d3==Contrast(iCont));
% 
%                 if isempty(TempIdx)==0
%                     d=ClmPLV_Trial(:,TempIdx);
%                     TempCorr=[];
%                     k=1;
%                     for iTrial=1:length(TempIdx)
%                         for jTrial=1:length(TempIdx)
%                             if iTrial<jTrial
%                                 TempCorr(k)=corr(d(:,iTrial),d(:,jTrial));
%                                 k=k+1;
%                             end
%                         end
%                     end
%                     AveCorrBetTrial(iResult,iCont,iSurr)=mean(TempCorr);
%                 else
%                     AveCorrBetTrial(iResult,iCont,iSurr)=NaN;
%                 end
%             end
%         end
% 
%     end
    
    SumData(iAnimal).AveCorrBetTrial=AveCorrBetTrial;
    
end
    
SumAveCorrBetTrial=[];
for iSurr=1:NSurr+1
    for iAnimal=1:NAnimal
        TempAveCorrBetTrial(:,:,iAnimal)=SumData(iAnimal).AveCorrBetTrial(:,:,iSurr);
    end
    SumAveCorrBetTrial(:,:,iSurr)=nanmean(TempAveCorrBetTrial,3);
end
    
% for iResult=1:4
%     for iCont=1:4
%         d=squeeze(SumAveCorrBetTrial(iResult,iCont,:));
%         if isnan(d(1))==0
%             p(iResult,iCont)=length(find(d(2:end)>=d(1)))/(NSurr-sum(isnan(d(2:end))));
%         else
%             p(iResult,iCont)=NaN;
%         end
%     end
% end
% p

p2=[];
for iResult=1:4
    for iCont=1:2
        d=squeeze(nanmean(SumAveCorrBetTrial(iResult,2*iCont-1:2*iCont,:),2));
        if isnan(d(1))==0
            p2(iResult,iCont)=length(find(d(2:end)>=d(1)))/(NSurr-sum(isnan(d(2:end))));
        else
            p2(iResult,iCont)=NaN;
        end
      
    end
end
p2
    
    
%ホルムーボンフェローニ補正後のp値(byChatGPT)
CorrPval=[0.728,0.000;0.104,0.584;0.128,0.072;0.104,0.000];

    
    
    
figure;
for iResult=1:4
    for iCont=1:2
        d=squeeze(nanmean(SumAveCorrBetTrial(iResult,2*iCont-1:2*iCont,:),2));
        subplot(2,4,iResult+4*(iCont-1));
        hist(d(2:end))
        h=findobj(gca,'Type','patch');
        set(h,'FaceColor',0.8*[1 1 1]);
        text(d(1),0,'↓','color',[0 0 1],'FontSize',30,'HorizontalAlignment','center','VerticalAlignment','baseline');
        xlim([min(d)-0.01,max(d)+0.01])
        xlabel('Corr','FontSize',12)
        title(['Result=',num2str(iResult),', Contrast=',num2str(iCont)],'FontSize',15)
    end
end

figure;
for iResult=1:4
    for iCont=1:4
        d=squeeze(nanmean(SumAveCorrBetTrial(iResult,iCont,:),2));
        subplot(4,4,(iCont-1)*4+iResult);
        hist(d(2:end))
        h=findobj(gca,'Type','patch');
        set(h,'FaceColor',0.8*[1 1 1]);
        text(d(1),0,'↓','color',[0 0 1],'FontSize',30,'HorizontalAlignment','center','VerticalAlignment','baseline');
        %xlim([min(d)-0.01,max(d)+0.01])
        xlabel('Corr','FontSize',12)
        title(['Result=',num2str(iResult),', Contrast=',num2str(iCont)],'FontSize',15)
    end
end

p3=[];
for iResult=1:4
    for iCont=1:4
        d=squeeze(nanmean(SumAveCorrBetTrial(iResult,iCont,:),2));
        if isnan(d(1))==0
            p3(iResult,iCont)=length(find(d(2:end)>=d(1)))/(NSurr-sum(isnan(d(2:end))));
        else
            p3(iResult,iCont)=NaN;
        end
    end
end

%% random labeling as the surrogate
% 分析二：通过打乱试验标签 (Shuffling Trial Labels) 进行统计检验
% 这种方法检验的是PLV模式与特定行为标签之间的关联是否显著。
% 它会为每个动物生成一个名为 AveCorrBetTrial_Shuffled 的新结果。

fprintf('\n--- 开始执行分析二：打乱试验标签 ---\n');

for iAnimal=1:NAnimal
    
    fprintf('正在处理动物 #%d (打乱标签分析)...\n', iAnimal);
    
    % 从 SumData 中获取所需数据
    ClmPLV_Trial_current = SumData(iAnimal).ClmPLV_Trial;
    ResultSeq_orig = SumData(iAnimal).ResultSeq;
    ContrastSeq_orig = SumData(iAnimal).ContrastSeq;
    
    % 如果没有柱内对，则跳过此动物
    if isempty(ClmPLV_Trial_current)
        fprintf('  动物 #%d 没有柱内对，跳过。\n', iAnimal);
        SumData(iAnimal).AveCorrBetTrial_Shuffled = nan(4, 4, NSurr + 1);
        continue;
    end
    
    AveCorrBetTrial_Shuffled = nan(4, 4, NSurr + 1); % 使用 NaN 初始化

    for iSurr=1:NSurr+1
        
        % iSurr=1 使用真实标签, iSurr>1 使用打乱后的标签
        if iSurr==1
            d2 = ResultSeq_orig;
            d3 = ContrastSeq_orig;
        else
            % 随机打乱试验结果和对比度的标签
            shuffled_indices = randperm(length(ResultSeq_orig));
            d2 = ResultSeq_orig(shuffled_indices);
            d3 = ContrastSeq_orig(shuffled_indices);
        end

        for iResult=1:4
            for iCont=1:4
                TempIdx = find(d2==iResult & d3==Contrast(iCont));

                % 只有当该分组的试验数大于1时，才能计算相关性
                if length(TempIdx) > 1
                    d = ClmPLV_Trial_current(:,TempIdx);
                    
                    % 使用 MATLAB 内置的 corr 函数直接计算相关性矩阵
                    corr_matrix = corr(d);
                    
                    % 提取上三角部分（不含对角线）并计算平均值
                    % tril(true(size(corr_matrix)),-1) 创建一个逻辑掩码
                    TempCorr = corr_matrix(tril(true(size(corr_matrix)),-1));
                    
                    if ~isempty(TempCorr)
                        AveCorrBetTrial_Shuffled(iResult,iCont,iSurr) = mean(TempCorr);
                    end
                end
            end
        end
    end
    
    SumData(iAnimal).AveCorrBetTrial_Shuffled = AveCorrBetTrial_Shuffled;
end

% --- 整合并可视化打乱标签后的结果 ---

% 整合所有动物的结果
SumAveCorrBetTrial_Shuffled=[];
for iSurr=1:NSurr+1
    TempAveCorrBetTrial_Shuffled = [];
    for iAnimal=1:NAnimal
        TempAveCorrBetTrial_Shuffled(:,:,iAnimal) = SumData(iAnimal).AveCorrBetTrial_Shuffled(:,:,iSurr);
    end
    SumAveCorrBetTrial_Shuffled(:,:,iSurr) = nanmean(TempAveCorrBetTrial_Shuffled, 3);
end

% 计算p值
p_shuffled = [];
for iResult=1:4
    for iCont=1:4
        d = squeeze(SumAveCorrBetTrial_Shuffled(iResult,iCont,:));
        if ~isnan(d(1))
            % 计算真实值在打乱分布中的排位
            p_shuffled(iResult,iCont) = sum(d(2:end) >= d(1)) / sum(~isnan(d(2:end)));
        else
            p_shuffled(iResult,iCont) = NaN;
        end
    end
end

fprintf('\n打乱标签分析得到的p值矩阵:\n');
disp(p_shuffled);

% 可视化
figure('Name', 'Shuffled Label Analysis Results');
sgtitle('打乱标签分析结果 (真实值 vs. 打乱分布)', 'FontWeight', 'bold');
for iResult=1:4
    for iCont=1:4
        d = squeeze(SumAveCorrBetTrial_Shuffled(iResult,iCont,:));
        subplot(4,4,(iCont-1)*4+iResult);
        
        if ~isnan(d(1))
            hist(d(2:end));
            h=findobj(gca,'Type','patch');
            set(h,'FaceColor',0.8*[1 1 1], 'EdgeColor', 'none');
            hold on;
            % 使用 xline 绘制一条垂线来标记真实值，更清晰
            xline(d(1), 'r-', 'LineWidth', 2);
            hold off;
            title(sprintf('R%d, C%d, p=%.3f', iResult, iCont, p_shuffled(iResult,iCont)));
        else
            title(sprintf('R%d, C%d, N/A', iResult, iCont));
            axis off;
        end
        
        if iResult == 1
            ylabel(sprintf('Contr %d', iCont));
        end
        if iCont == 4
            xlabel(sprintf('Res %d', iResult));
        end
    end
end
%% 分析三: 将Contrast分组后再次进行打乱标签分析
% 这个部分利用上一步“打乱标签分析”的结果，将4个对比度合并为高/低两组，
% 然后重新计算p值并可视化，以检查是否存在更概括性的显著结果。

fprintf('\n--- 分析三: 将Contrast分为高/低两组后重新计算p值 ---\n');

% 初始化p值矩阵 (4个结果 x 2个对比度组)
p_shuffled_grouped = nan(4, 2);
% 用于存储分组后的数据以供绘图
SumAveCorrBetTrial_Shuffled_Grouped = nan(4, 2, NSurr + 1);

for iResult = 1:4
    for iContGroup = 1:2 % 1: 低对比度组, 2: 高对比度组
        
        % 定义每个组包含的原始对比度索引
        % iContGroup=1 -> 对应原始的 Contrast 1,2 (0.001, 0.01)
        % iContGroup=2 -> 对应原始的 Contrast 3,4 (0.1, 1)
        contrast_indices = (2*iContGroup-1):(2*iContGroup);
        
        % 从已有的结果中提取数据，并沿第2维（对比度维度）求平均
        d_grouped = squeeze(nanmean(SumAveCorrBetTrial_Shuffled(iResult, contrast_indices, :), 2));
        
        % 存储分组后的数据，方便后续绘图
        SumAveCorrBetTrial_Shuffled_Grouped(iResult, iContGroup, :) = d_grouped;
        
        % 计算p值：比较真实值(d_grouped(1))与打乱后的分布(d_grouped(2:end))
        if ~isnan(d_grouped(1))
            p_shuffled_grouped(iResult, iContGroup) = sum(d_grouped(2:end) >= d_grouped(1)) / sum(~isnan(d_grouped(2:end)));
        else
            p_shuffled_grouped(iResult, iContGroup) = NaN;
        end
    end
end

fprintf('\n打乱标签分析 (Contrast分组后) 得到的p值矩阵:\n');
disp(p_shuffled_grouped);

% --- 可视化分组后的结果 ---
figure('Name', 'Shuffled Label Analysis Results (Grouped Contrast)');
sgtitle('打乱标签分析结果 (Contrast分组)', 'FontWeight', 'bold');
ContrastLabels = {'Low Contrast (<10%)', 'High Contrast (>=10%)'};

for iResult = 1:4
    for iContGroup = 1:2
        d = squeeze(SumAveCorrBetTrial_Shuffled_Grouped(iResult, iContGroup, :));
        subplot(2, 4, iResult + 4*(iContGroup-1));
        
        if ~isnan(d(1))
            hist(d(2:end));
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none');
            hold on;
            xline(d(1), 'r-', 'LineWidth', 2);
            hold off;
            title(sprintf('Res %d, p=%.4f', iResult, p_shuffled_grouped(iResult, iContGroup)));
            xlim_vals = xlim;
            ylim_vals = ylim;
            text(xlim_vals(1), ylim_vals(2)*1.1, ContrastLabels{iContGroup}, 'VerticalAlignment', 'bottom');
        else
            title(sprintf('Res %d, N/A', iResult));
            axis off;
        end
    end
end
%% 解码分析 (改进版): 解码感知结果与行为决策
% 目标: 在高对比度试验中，使用PLV特征分别解码：
% 1. 感知结果 (Correct vs. Incorrect)
% 2. 行为决策 (Lick vs. No Lick)

fprintf('\n--- 开始执行解码分析 (改进版) ---\n');

% --- 1. 整合所有动物的PLV特征数据 (通过PCA降维) ---
% 这部分与之前相同，汇集所有数据以便后续筛选
N_COMPONENTS = 100;
All_X_pca = [];
All_ResultSeq = [];
All_ContrastSeq = [];

for iAnimal = 1:NAnimal
    animal_name = erase(SumData(iAnimal).OriginalFileName, {'ROISegTraceTable_', '.mat'});
    fprintf('正在处理动物: %s\n', animal_name);
    
    if isempty(SumData(iAnimal).ClmPLV_Trial)
        fprintf('  警告: 动物 %s 没有找到 in-column pairs，跳过此动物。\n', animal_name);
        continue;
    end
    
    X_animal = SumData(iAnimal).ClmPLV_Trial';
    
    num_features = size(X_animal, 2);
    if num_features < N_COMPONENTS
        current_n_components = num_features;
    else
        current_n_components = N_COMPONENTS;
    end
    
    [~, score] = pca(X_animal, 'NumComponents', current_n_components);
    
    if current_n_components < N_COMPONENTS
        X_animal_pca = zeros(size(score, 1), N_COMPONENTS);
        X_animal_pca(:, 1:current_n_components) = score;
    else
        X_animal_pca = score;
    end

    All_X_pca = [All_X_pca; X_animal_pca];
    All_ResultSeq = [All_ResultSeq; SumData(iAnimal).ResultSeq];
    All_ContrastSeq = [All_ContrastSeq; SumData(iAnimal).ContrastSeq];
end

fprintf('所有动物数据整合完成。\n');

% --- 2. 筛选高对比度试验并创建标签 ---
high_contrast_idx = find((All_ContrastSeq == 0.1) | (All_ContrastSeq == 1));

% 仅使用高对比度试验的数据
X_high_contrast = All_X_pca(high_contrast_idx, :);
Result_high_contrast = All_ResultSeq(high_contrast_idx);

% 创建两个解码任务的标签
% 任务一：感知结果 (Correct vs. Incorrect)
Y_perception = (Result_high_contrast == 1) | (Result_high_contrast == 4); % 1=Correct, 0=Incorrect

% 任务二：行为决策 (Lick vs. No Lick)
Y_choice = (Result_high_contrast == 1) | (Result_high_contrast == 3); % 1=Lick, 0=No Lick

% --- 3. 训练并评估“感知结果”解码器 ---
fprintf('\n--- 任务一: 解码感知结果 (Correct vs. Incorrect) ---\n');
fprintf('  Correct 试验数: %d\n', sum(Y_perception==1));
fprintf('  Incorrect 试验数: %d\n', sum(Y_perception==0));

CV_FOLDS = 5;
t = templateSVM('KernelFunction', 'linear', 'Standardize', true, 'BoxConstraint', 0.1);
CVMdl_perception = fitcecoc(X_high_contrast, Y_perception, 'Learners', t, 'Crossval', 'on', 'KFold', CV_FOLDS);
accuracy_perception = 1 - kfoldLoss(CVMdl_perception, 'LossFun', 'ClassifError');
fprintf('  交叉验证平均准确率: %.3f\n', accuracy_perception);

% 可视化
figure('Name', 'Perception Decoder (Correct vs. Incorrect)');
predictions_perception = kfoldPredict(CVMdl_perception);
cm_perception = confusionchart(Y_perception, predictions_perception, ...
    'Title', '解码器一: 感知结果', ...
    'RowSummary', 'row-normalized');
%cm_perception.ClassLabels = {'Incorrect', 'Correct'};


% --- 4. 训练并评估“行为决策”解码器 ---
fprintf('\n--- 任务二: 解码行为决策 (Lick vs. No Lick) ---\n');
fprintf('  Lick 试验数: %d\n', sum(Y_choice==1));
fprintf('  No Lick 试验数: %d\n', sum(Y_choice==0));

CVMdl_choice = fitcecoc(X_high_contrast, Y_choice, 'Learners', t, 'Crossval', 'on', 'KFold', CV_FOLDS);
accuracy_choice = 1 - kfoldLoss(CVMdl_choice, 'LossFun', 'ClassifError');
fprintf('  交叉验证平均准确率: %.3f\n', accuracy_choice);

% 可视化
figure('Name', 'Choice Decoder (Lick vs. No Lick)');
predictions_choice = kfoldPredict(CVMdl_choice);
cm_choice = confusionchart(Y_choice, predictions_choice, ...
    'Title', '解码器二: 行为决策', ...
    'RowSummary', 'row-normalized');
%cm_choice.ClassLabels = {'No Lick', 'Lick'};

disp('基于PLV特征的双解码任务分析完成。');

%% 解码分析三: 解码准确率随神经元对数量的变化
% 目标: 观察随着纳入分析的in-column pair数量增加，
%       “感知结果”和“行为决策”解码器的准确率如何变化。
%
% 1. 整合所有动物的原始PLV数据（不使用PCA）。
% 2. 筛选出高对比度试验。
% 3. 在一个循环中，逐步增加使用的神经元对数量。
% 4. 在每一步，分别训练和评估两个解码器，并记录准确率。
% 5. 绘制准确率曲线图。

fprintf('\n--- 开始执行解码分析三: 准确率 vs. 神经元对数量 ---\n');

% --- 1. 整合所有动物的原始PLV数据 ---
all_animals_plv = [];
all_animals_results = [];
all_animals_contrasts = [];

for iAnimal = 1:NAnimal
    if isempty(SumData(iAnimal).ClmPLV_Trial)
        continue;
    end
    
    % 获取PLV数据 (转置为 试验数 x 神经元对数)
    plv_data = SumData(iAnimal).ClmPLV_Trial';
    
    % 获取行为和对比度数据
    result_seq = SumData(iAnimal).ResultSeq;
    contrast_seq = SumData(iAnimal).ContrastSeq;
    
    % 确保数据行数匹配
    if size(plv_data, 1) ~= length(result_seq)
        fprintf('  警告: 动物 #%d 的PLV数据与行为数据行数不匹配，跳过。\n', iAnimal);
        continue;
    end
    
    % 将当前动物的数据添加到总数据池中
    % 注意：这里我们简单地将所有动物的神经元对横向拼接
    % 这假设了每个神经元对都是一个独立的特征来源
    all_animals_plv = blkdiag(all_animals_plv, plv_data);
    all_animals_results = [all_animals_results; result_seq];
    all_animals_contrasts = [all_animals_contrasts; contrast_seq];
end

% --- 2. 筛选高对比度试验 ---
high_contrast_idx = find((all_animals_contrasts == 0.1) | (all_animals_contrasts == 1));
X_full = all_animals_plv(high_contrast_idx, :);
Result_full = all_animals_results(high_contrast_idx);

% 创建两个解码任务的标签
Y_perception = (Result_full == 1) | (Result_full == 4); % 1=Correct, 0=Incorrect
Y_choice = (Result_full == 1) | (Result_full == 3);     % 1=Lick, 0=No Lick

% --- 3. 循环增加神经元对数量并进行解码 ---
num_total_pairs = size(X_full, 2);
STEP_SIZE = 5; % 每次增加5个神经元对
pair_counts = 1:STEP_SIZE:num_total_pairs;
if pair_counts(end) ~= num_total_pairs % 确保包含所有神经元对的情况
    pair_counts(end+1) = num_total_pairs;
end

accuracy_curve_perception = zeros(size(pair_counts));
accuracy_curve_choice = zeros(size(pair_counts));

fprintf('开始循环解码，总共 %d 个神经元对，步长为 %d...\n', num_total_pairs, STEP_SIZE);
tic;
parfor i = 1:length(pair_counts)
    num_pairs = pair_counts(i);
    fprintf('  正在处理: 使用 %d 个神经元对...\n', num_pairs);
    
    % 随机选择指定数量的神经元对作为特征
    % 为了结果稳定，每次循环都应随机选择，但为了演示，这里我们按顺序选择
    % 如果要获得更稳健的结果，应多次重复随机选择并取平均
    feature_indices = 1:num_pairs;
    X_subset = X_full(:, feature_indices);
    
    % 定义SVM模板
    t = templateSVM('KernelFunction', 'linear', 'Standardize', true, 'BoxConstraint', 0.1);
    CV_FOLDS = 5;
    
    % 训练感知解码器
    CVMdl_p = fitcecoc(X_subset, Y_perception, 'Learners', t, 'Crossval', 'on', 'KFold', CV_FOLDS);
    accuracy_curve_perception(i) = 1 - kfoldLoss(CVMdl_p, 'LossFun', 'ClassifError');
    
    % 训练行为决策解码器
    CVMdl_c = fitcecoc(X_subset, Y_choice, 'Learners', t, 'Crossval', 'on', 'KFold', CV_FOLDS);
    accuracy_curve_choice(i) = 1 - kfoldLoss(CVMdl_c, 'LossFun', 'ClassifError');
end
toc;

% --- 4. 可视化结果 ---
figure('Name', 'Decoding Accuracy vs. Number of Pairs');
plot(pair_counts, accuracy_curve_perception, 'b-o', 'LineWidth', 2, 'DisplayName', '感知结果 (Correct/Incorrect)');
hold on;
plot(pair_counts, accuracy_curve_choice, 'r-s', 'LineWidth', 2, 'DisplayName', '行为决策 (Lick/No Lick)');
hold off;

% 计算机会水平
chance_perception = max(mean(Y_perception), 1-mean(Y_perception));
chance_choice = max(mean(Y_choice), 1-mean(Y_choice));
yline(chance_perception, 'b--', 'DisplayName', sprintf('感知机会水平 (%.2f)', chance_perception));
yline(chance_choice, 'r--', 'DisplayName', sprintf('决策机会水平 (%.2f)', chance_choice));

title('解码准确率随神经元对数量的变化');
xlabel('使用的In-Column Pair数量');
ylabel('交叉验证准确率');
legend('show', 'Location', 'southeast');
grid on;
ylim([0 1]); % 准确率范围在0到1之间

disp('解码准确率曲线分析完成。');

%% 编码模型 (改进版): 基于行为结果定义 Target/Non-Target
% 目标: 检验呈现"目标刺激"(导致Hit/Miss)与呈现"非目标刺激"(导致FA/CR)
%       是否引起了不同的PLV响应强度。
%
% 1. 筛选出所有高对比度(>=10%)的试验。
% 2. 在这些试验中，根据行为结果分为"Target"和"Non-Target"两组。
% 3. 对每组，计算其所有试验的平均柱内对PLV。
% 4. 使用t-test比较两组的平均PLV是否有显著差异。

fprintf('\n--- 开始执行编码模型分析 (基于行为定义) ---\n');

% 初始化用于存储所有动物数据的变量
all_mean_plv = [];
all_group_labels = []; % 1 for Target, 2 for Non-Target

for iAnimal = 1:NAnimal
    animal_name = erase(SumData(iAnimal).OriginalFileName, {'ROISegTraceTable_', '.mat'});
    
    % 检查是否存在柱内对PLV数据
    if isempty(SumData(iAnimal).ClmPLV_Trial)
        fprintf('  动物 %s 没有找到 in-column pairs，跳过此动物的编码模型分析。\n', animal_name);
        continue;
    end
    
    % 计算每个试验的平均PLV
    mean_plv_per_trial = mean(SumData(iAnimal).ClmPLV_Trial, 1);
    
    % 获取行为和对比度序列
    result_seq = SumData(iAnimal).ResultSeq;
    contrast_seq = SumData(iAnimal).ContrastSeq;
    
    % 1. 创建高对比度试验的筛选掩码
    high_contrast_mask = (contrast_seq == 0.1) | (contrast_seq == 1);
    
    % 2. 创建行为学定义的 Target/Non-Target 掩码
    target_trial_mask = (result_seq == 1) | (result_seq == 2); % Hit or Miss
    non_target_trial_mask = (result_seq == 3) | (result_seq == 4); % FA or CR
    
    % 结合两个掩码，找到最终要分析的试验索引
    target_indices = find(target_trial_mask & high_contrast_mask);
    non_target_indices = find(non_target_trial_mask & high_contrast_mask);
    
    % 提取对应分组的平均PLV值
    target_plv = mean_plv_per_trial(target_indices);
    non_target_plv = mean_plv_per_trial(non_target_indices);
    
    % 汇集所有动物的数据
    all_mean_plv = [all_mean_plv, target_plv, non_target_plv];
    all_group_labels = [all_group_labels, ones(1, length(target_plv)), 2*ones(1, length(non_target_plv))];
end

% 3. 对整合后的所有数据进行统计检验
fprintf('正在对所有高对比度试验的整合数据进行t-test...\n');
target_data = all_mean_plv(all_group_labels == 1);
non_target_data = all_mean_plv(all_group_labels == 2);

[h, p_value, ci, stats] = ttest2(target_data, non_target_data);

fprintf('t-test 结果:\n');
fprintf('  p-value: %.5f\n', p_value);
if h
    fprintf('  结论: Target组和Non-Target组的平均PLV存在显著差异。\n');
else
    fprintf('  结论: 未发现两组平均PLV存在显著差异。\n');
end
fprintf('  t-statistic: %.3f\n', stats.tstat);
fprintf('  Target 组 (Hit/Miss) 平均PLV: %.4f (n=%d trials)\n', mean(target_data), length(target_data));
fprintf('  Non-Target 组 (FA/CR) 平均PLV: %.4f (n=%d trials)\n', mean(non_target_data), length(non_target_data));


% 4. 可视化结果
figure('Name', 'Encoding Model (Behaviorally Defined)');
boxplot(all_mean_plv, all_group_labels, 'Labels', {'Target Trials (Hit/Miss)', 'Non-Target Trials (FA/CR)'});
title({'编码模型: PLV响应强度', '(仅包含高对比度试验)'});
ylabel('平均柱内对PLV (Mean In-Column PLV)');
grid on;

disp('基于行为定义的编码模型分析完成。');