ClmDis=10;
AdClmDis=30;
ExDis=20;
Th=1;
StartBin=1;
EndBin=20;
Contrast=[0.001,0.01,0.1,1];
NTrial=size(SumData(1).Data.Segmented_trace{1,1},1);
NAnimal=3;
NSurr=1000;

NTimePoints=EndBin;
% 時間軸
Fs = 2.37; % サンプリング周波数 [Hz]
time = (0:NTimePoints-1) / Fs; % 時間 [s]


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