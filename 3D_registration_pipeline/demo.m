%E14?解析用 - 高速化版
EdgeUM=30;
ZlimMax=3;
NData=3;
RangeUM=[7,10]; % tangential距離フィルタ
WithinUM=30;

%出力データの保存先
disp('==保存先フォルダを選択してください==')
OutDir=uigetdir;
cd(OutDir)

%Surrogate dataの作成
NSurr=0;

for iData=1:NData
    for iSurr=1:NSurr
        EdUPosUM=Data.EdU{iData,1}.PosUM;
        NEdUCell=size(EdUPosUM,1);
        FTPosUM=Data.FlashTag{iData,1}.PosUM;
        AllPosUM=[EdUPosUM;FTPosUM];
        NAllCell=size(AllPosUM,1);

        TempEduIdx=randperm(NAllCell,NEdUCell);
        TempFTIdx = setdiff(1:NAllCell, TempEduIdx);

        Data.EdU{iData,iSurr+1}.PosUM=AllPosUM(TempEduIdx,:);
        Data.FT{iData,iSurr+1}.PosUM=AllPosUM(TempFTIdx,:);
    end
end

tic
for iSurr=1:NSurr+1 %数字注意

%まず細胞座標を近似平面がxy平面となるように回転させる。
for iData=1:NData

    PosUM0=Data.EdU{iData,iSurr}.PosUM;  %元の細胞座標
    ImgSizeUM=Data.EdU{iData,1}.ImgSizeUM;
    d=PosUM0;
    f = fit([d(:,1), d(:,2)],d(:,3),'poly11');
    
    % 図の作成をスキップ（高速化）
    % figure;
    % plot(f,[d(:,1), d(:,2)],d(:,3))
    % axis equal
    % close all

    x=0:ceil(max(d(:,1)));
    y=0:ceil(max(d(:,2)));
    [X,Y]=meshgrid(x,y);
    Z=f(X,Y);

    NVec=[f.p10,f.p01,-1];

    %%%
    ClmAxisVec=NVec;
    ClmAxisVec(3)=sqrt(1-(ClmAxisVec(1)^2+ClmAxisVec(2)^2));
    [tgtrng,tgtang] = rangeangle(ClmAxisVec');

    TempAngle=[tgtang(1),tgtang(2)-90];

    Theta2=deg2rad(-1*TempAngle(1));
    Phi2=deg2rad(-1*TempAngle(2));
    ClmPrjMtx1=[cos(Theta2) -sin(Theta2);sin(Theta2) cos(Theta2)];
    ClmPrjMtx2=[cos(Phi2) -sin(Phi2);sin(Phi2) cos(Phi2)];

    TempVec=ClmAxisVec;

    Temp4=ClmPrjMtx1*[TempVec(1);TempVec(2)];
    TempVec(1)=Temp4(1);
    TempVec(2)=Temp4(2);
    Temp3=ClmPrjMtx2*[TempVec(1);TempVec(3)];
    TempVec(1)=Temp3(1);
    TempVec(3)=Temp3(2);
    [temprng,tempang] = rangeangle(TempVec');

    disp('ClmAxisVecの回転後の行列')
    TempVec  %法線ベクトルがz軸と平行になったかを確認する。

    %%%
    RefDim1Range=[EdgeUM,ImgSizeUM(2)-EdgeUM];
    RefDim2Range=[EdgeUM,ImgSizeUM(1)-EdgeUM];

    RefPos=find(PosUM0(:,1)>=RefDim2Range(1)&PosUM0(:,1)<=RefDim2Range(2)&PosUM0(:,2)>=RefDim1Range(1)&PosUM0(:,2)<=RefDim1Range(2));
    TempPosUM1=PosUM0(RefPos,:);
    TempPosUM2=PosUM0;
    RefPos0{1,iData}=RefPos; %RefPosを保存

    % 【最適化1】ベクトル化された座標回転
    tTempPosUM1 = applyRotationVectorized(TempPosUM1, ClmPrjMtx1, ClmPrjMtx2);
    tTempPosUM2 = applyRotationVectorized(TempPosUM2, ClmPrjMtx1, ClmPrjMtx2);

    FlatEdUPosUM{1,iData}=tTempPosUM2;
end

%上で得られたFlatPosUMを元に、さまざまな方向に振って解析をおこなう。

% % 結果を保存する一時ディレクトリ
% tempDir = '/Users/marubook2/Documents/MATLAB';
% cd(tempDir);

Angle1=0:3:359;
Angle2=-30:1:30;

% 【最適化2】単位ベクトルの事前計算
UnitVec=[0 0 1];
UnitVec3 = precomputeUnitVectors(Angle1, Angle2, UnitVec);

% 【最適化3】並列化可能な構造に変更
nTotal=size(UnitVec3,1)*NData;
AllExrPosUM = cell(size(UnitVec3,1), NData);

h1 = waitbar(0, '回転座標計算中...'); % waitbar作成
for iData=1:NData
    % 【最適化4】parforループ（並列計算ツールボックスがある場合）
    for iRot=1:size(UnitVec3,1)
        ClmAxisVec=UnitVec3(iRot,:);
        ClmAxisVec(3)=sqrt(1-(ClmAxisVec(1)^2+ClmAxisVec(2)^2));
        [tgtrng,tgtang] = rangeangle(ClmAxisVec');

        TempAngle=[tgtang(1),tgtang(2)-90];

        Theta2=deg2rad(-1*TempAngle(1));
        Phi2=deg2rad(-1*TempAngle(2));
        ClmPrjMtx1=[cos(Theta2) -sin(Theta2);sin(Theta2) cos(Theta2)];
        ClmPrjMtx2=[cos(Phi2) -sin(Phi2);sin(Phi2) cos(Phi2)];

        TempPosUM1=FlatEdUPosUM{1,iData}(RefPos0{1,iData},:);
        TempPosUM2=FlatEdUPosUM{1,iData};

        % ベクトル化された回転
        tTempPosUM1 = applyRotationVectorized(TempPosUM1, ClmPrjMtx1, ClmPrjMtx2);
        tTempPosUM2 = applyRotationVectorized(TempPosUM2, ClmPrjMtx1, ClmPrjMtx2);

        % 【最適化5】相対座標計算の高速化
        ExrPosUM = computeRelativePositionsVectorized(tTempPosUM1, tTempPosUM2, ZlimMax);
        TempIdx=find(abs(ExrPosUM(:,1))<=WithinUM&abs(ExrPosUM(:,2))<=WithinUM&abs(ExrPosUM(:,3))<=ZlimMax);
         AllExrPosUM{iRot,iData} = ExrPosUM(TempIdx,:);
%         IndExrPosUM{iRot,1} = ExrPosUM;
        
        clear tTempPosUM1 tTempPosUM2 ExrPosUM TempIdx
        
        if mod((iData-1)*size(UnitVec3,1) + iRot, 10)==0
            waitbar(((iData-1)*size(UnitVec3,1) + iRot) / nTotal, h1, ...
            sprintf('回転座標計算中 (%d/%d)', (iData-1)*size(UnitVec3,1) + iRot, nTotal));
            toc
        end
        
    end
    
%     % 保存
%     saveFile = fullfile(tempDir, sprintf('IndExrPosUM_%d.mat', iData));
%     save(saveFile, 'IndExrPosUM', '-v7.3');
%     clear IndExrPosUM
%     
end
close(h1)

% for iData=1:NData
%     load(['IndExrPosUM_',num2str(iData),'.mat']);
%     for iRot=1:size(IndExrPosUM,1)
%         AllExrPosUM{iRot,iData}=IndExrPosUM{iRot,1};
%     end
%     clear IndExrPosUM
%     tocf(2)
% end

disp('回転座標計算終了')

% 【最適化6】ガウシアンフィルタの事前計算
Sigma=4;
BinWidth=0.3;
Bin1=-30:BinWidth:30;
Bin2=-30:BinWidth:30;

% ガウシアンカーネルの事前計算
gaussKernel = precomputeGaussianKernel(Sigma, BinWidth);

nTotal=size(UnitVec3,1)*NData;
h1 = waitbar(0, '細胞密度計算中...'); % waitbar作成
TempColorMap=zeros(length(Bin1), length(Bin2), NData, size(UnitVec3,1));

for iData=1:NData
    for iRot=1:size(UnitVec3,1)
        % 【最適化7】高速化されたガウシアンフィルタ
%         GaussMap = GaussConvOptimized(AllExrPosUM{iRot,iData}(:,1:2), Bin1, Bin2, gaussKernel);
%         load(['IndExrPosUM_',num2str(iData),'.mat']);
%         GaussMap = GaussConv(IndExrPosUM{iRot,1}(:,1:2),Bin1,Bin2,Sigma);
%         TempColorMap(:,:,iData,iRot)=GaussMap;
        
        GaussMap = GaussConv(AllExrPosUM{iRot,iData}(:,1:2),Bin1,Bin2,Sigma);
        TempColorMap(:,:,iData,iRot)=GaussMap;

        if mod((iData-1)*size(UnitVec3,1) + iRot, 10)==0
            waitbar(((iData-1)*size(UnitVec3,1) + iRot) / nTotal, h1, ...
                sprintf('細胞密度計算中 (%d/%d)', (iData-1)*size(UnitVec3,1) + iRot, nTotal));
        end
        toc
    end
end
close(h1)
disp('細胞密度計算終了')

% clear AllExrPosUM

% 【最適化8】極座標変換の事前計算
[polarGrid, rangeMask] = precomputePolarGrid(Bin1, Bin2, RangeUM, BinWidth);

%tangential距離フィルタ内の細胞密度を極座標の波形として算出する。
NAngle = 360;
OutData=cell(1,NData);
SumrCellDensity=nan(NAngle,NData,size(UnitVec3,1));

h1 = waitbar(0, 'Tangential範囲内密度計算中...'); % waitbar作成
nTotal = NData * size(UnitVec3,1);

for iData = 1:NData
    for iRot = 1:size(UnitVec3,1)
        % 【最適化9】ベクトル化された極座標密度計算
        rCellDensity = computePolarDensityVectorized(TempColorMap(:,:,iData,iRot), polarGrid, rangeMask, NAngle);
        SumrCellDensity(:,iData,iRot) = rCellDensity;
        
        % 【最適化10】高速化されたフーリエ変換
%         power_spec = computeCircularPowerSpectrum(rCellDensity, 10);

        power_spec=[];
        theta_deg2 = 0:359;
        theta_rad = 2 * pi * theta_deg2 / 360;
        max_harmonic = 10;  % 1-10回の対称性まで調べる

        y = SumrCellDensity(:,iData,iRot)';
        N = length(y);

        for k = 1:max_harmonic
            A_n = sum(y .* exp(-1i * k * theta_rad)) / N;
            power_spec(k) = abs(A_n).^2;
        end
        OutData{1,iData}(iRot,:) = power_spec;
        
        if mod((iData-1)*size(UnitVec3,1) + iRot, 10)==0
            waitbar(((iData-1)*size(UnitVec3,1) + iRot) / nTotal, h1, ...
            sprintf('密度計算中 (%d/%d)', (iData-1)*size(UnitVec3,1) + iRot, nTotal));
            toc
        end
    end
end
close(h1);
disp('Tangential範囲内密度計算終了')

% 残りの処理は元のコードと同じ...
HexaPower=zeros(NData,size(UnitVec3,1));
% for iData=1:NData
%     for iRot=1:size(UnitVec3,1)
%         d=OutData{1,iData}(iRot,:);
%         HexaPower(iData,iRot)=d(6)/sum(d);
%     end
% end


for iData=1:NData
    for iRot=1:size(UnitVec3,1)
        d=TempColorMap(:,:,iData,iRot);
        LocalMax = imregionalmax(d);
        NPeaks(iData,iRot)=length(find((LocalMax.*rangeMask)==1));
        if NPeaks(iData,iRot)==6
            d=OutData{1,iData}(iRot,:);
            HexaPower(iData,iRot)=(d(6)/sum(d));
%             HexaPower(iData,iRot)=std(SumrCellDensity(:,iData,iRot));
        else
            HexaPower(iData,iRot)=0;
        end
    end
end


Top1Vec = zeros(NData,3);
for iData=1:NData
    threshold = prctile(HexaPower(iData,:), 99.9);
    top1_idx = find(HexaPower(iData,:) >= threshold);
    Top1Vec(iData,:)=mean(UnitVec3(top1_idx,:),1);
end

% 可視化
figure
for iData=1:NData
    subplot(1,NData,iData)
    hold on;
    scatter(UnitVec3(:,1), UnitVec3(:,2), 50, HexaPower(iData,:), 'filled');
    plot(Top1Vec(iData,1),Top1Vec(iData,2),'r+','MarkerSize',20)
    title('Hexagon Power', 'FontSize', 12);
    colorbar;
    axis equal
end

% 最終処理
FinalColorMap=zeros(length(Bin1), length(Bin2), NData);
FinalPosUM=cell(1,NData);
FinalRelatPosUM=cell(1,NData);

for iData=1:NData
    ClmAxisVec=Top1Vec(iData,:);
    ClmAxisVec(3)=sqrt(1-(ClmAxisVec(1)^2+ClmAxisVec(2)^2));
    [tgtrng,tgtang] = rangeangle(ClmAxisVec');

    TempAngle=[tgtang(1),tgtang(2)-90];

    Theta2=deg2rad(-1*TempAngle(1));
    Phi2=deg2rad(-1*TempAngle(2));
    ClmPrjMtx1=[cos(Theta2) -sin(Theta2);sin(Theta2) cos(Theta2)];
    ClmPrjMtx2=[cos(Phi2) -sin(Phi2);sin(Phi2) cos(Phi2)];

    TempPosUM1=FlatEdUPosUM{1,iData}(RefPos0{1,iData},:);
    TempPosUM2=FlatEdUPosUM{1,iData};

    tTempPosUM1 = applyRotationVectorized(TempPosUM1, ClmPrjMtx1, ClmPrjMtx2);
    tTempPosUM2 = applyRotationVectorized(TempPosUM2, ClmPrjMtx1, ClmPrjMtx2);

    FinalPosUM{1,iData}=tTempPosUM2;
    
    ExrPosUM = computeRelativePositionsVectorized(tTempPosUM1, tTempPosUM2, ZlimMax);
    FinalRelatPosUM{1,iData}=ExrPosUM;
    
%     GaussMap = GaussConvOptimized(ExrPosUM(:,1:2), Bin1, Bin2, gaussKernel);
    GaussMap = GaussConv(AllExrPosUM{iRot,iData}(:,1:2),Bin1,Bin2,Sigma);
    FinalColorMap(:,:,iData)=GaussMap;
end

% 極座標表示の細胞密度波形を計算する
FinalrCellDensity = zeros(NAngle, NData);
FinalPower_Spec = cell(1,NData);

for iData = 1:NData
    rCellDensity = computePolarDensityVectorized(FinalColorMap(:,:,iData), polarGrid, rangeMask, NAngle);
    FinalrCellDensity(:,iData) = rCellDensity;
    
%     power_spec = computeCircularPowerSpectrum(rCellDensity, 10);

        %Circular power spectrum analysis

    power_spec=[];
    theta_deg2 = 0:359;
    theta_rad = 2 * pi * theta_deg2 / 360;
    max_harmonic = 10;  % 1-10回の対称性まで調べる

    y = SumrCellDensity(:,iData,iRot)';
    N = length(y);

    for k = 1:max_harmonic
        A_n = sum(y .* exp(-1i * k * theta_rad)) / N;
        power_spec(k) = abs(A_n).^2;
    end
    
    FinalPower_Spec{1,iData}=power_spec;
end

% 回転補正
% RotCorrect=zeros(1,NData);
% d1=FinalrCellDensity(:,1)/mean(FinalrCellDensity(:,1))-1;
% d2=FinalrCellDensity(:,2)/mean(FinalrCellDensity(:,2))-1;
% if NData >= 3
%     d3=FinalrCellDensity(:,3)/mean(FinalrCellDensity(:,3))-1;
% end
% 
% r=xcorr(smooth(d1),smooth(d2));
% [~,Idx1]=max(r);
% RotCorrect(1) = abs(Idx1-360);
% 
% if NData >= 3
%     r=xcorr(smooth(d3),smooth(d2));
%     [~,Idx2]=max(r);
%     RotCorrect(3) = abs(Idx2-360);
% end

iData=1;
img1=IndividualColorMap(:,:,iData);
[x1,y1]=find(imregionalmax(img1).*rangeMask);
NpeakPosi1=[x1,y1];
Npeaks1=length(find(imregionalmax(img1).*rangeMask));

RotCorrect(1)=0;

for jData=1:NData-1
    
    iData=jData+1;
    img2=IndividualColorMap(:,:,iData);

    TempSumDis=nan(360,1);
    for iAngle=1:360
        TempImg=imrotate(img2,iAngle-1,'crop');
        [x2,y2]=find(imregionalmax(TempImg).*rangeMask);
        NpeakPosi2=[x2,y2];
        Npeaks2=length(find(imregionalmax(TempImg).*rangeMask));

        TempTable=nan(Npeaks1,Npeaks2);
        for pk1=1:Npeaks1
            for pk2=1:Npeaks2
                TempTable(pk1,pk2)=norm(NpeakPosi1(pk1,:)-NpeakPosi2(pk2,:));
            end
        end

        mpk=min([Npeaks1,Npeaks2]);

        TempDisTable=nan(mpk,1);
        if Npeaks1>=Npeaks2
            for pk=1:mpk
                TempDisTable(pk)=min(TempTable(pk,:));
            end
        else
            for pk=1:mpk
                TempDisTable(pk)=min(TempTable(:,pk));
            end  
        end

        TempSumDis(iAngle)=sum(TempDisTable);

    end

    TempSumDis2=TempSumDis;
    RotCorrect(jData+1)=find(TempSumDis2==min(TempSumDis2))-1-360;

end


IndividualColorMap = zeros(length(Bin1), length(Bin2), NData);

for iData=1:NData
    TemprPosUM=FinalRelatPosUM{1,iData};
    
    % 回転補正の適用
%     CorrRotMtx=[cosd(RotCorrect(iData)) -sind(RotCorrect(iData));
%                 sind(RotCorrect(iData)) cosd(RotCorrect(iData))];
CorrRotMtx=[cos(ang2rad(RotCorrect(iData))) -sin(ang2rad(RotCorrect(iData)));sin(ang2rad(RotCorrect(iData))) cos(ang2rad(RotCorrect(iData)))];

    TemprPosUM2 = TemprPosUM;
    TemprPosUM2(:,1:2) = (CorrRotMtx * TemprPosUM(:,1:2)')';
    
%     IndividualColorMap(:,:,iData) = GaussConvOptimized(TemprPosUM2(:,1:2), Bin1, Bin2, gaussKernel);


    Temp=[];
    Temp(:,:,1)=GaussConv(TemprPosUM(:,1:2),Bin1,Bin2,Sigma);%回転前
    Temp(:,:,2)=GaussConv(TemprPosUM2(:,1:2),Bin1,Bin2,Sigma);%回転後

    figure;
    subplot(1,2,1);
    imagesc(Temp(:,:,1));
    caxis([mean(mean(Temp(:,:,1))),max(max(Temp(:,:,1)))]);
    axis equal
    subplot(1,2,2);
    imagesc(Temp(:,:,2));
    caxis([mean(mean(Temp(:,:,2))),max(max(Temp(:,:,2)))]);
    axis equal


    IndividualColorMap(:,:,iData)=Temp(:,:,2);

end

% 合計カラー画像の計算
SumColorMap = sum(IndividualColorMap,3);

figure
for iData=1:NData
    subplot(floor(sqrt(NData+1)),ceil(sqrt(NData+1)),iData);
    imagesc(IndividualColorMap(:,:,iData));
%     caxis([min(FinalrCellDensity(:,iData)),max(FinalrCellDensity(:,iData))])
    dd=IndividualColorMap(:,:,iData).*rangeMask;
    dd(dd==0)=nan;
    caxis([min(min(dd)),max(max(dd))])
    axis equal
    colorbar
end
subplot(floor(sqrt(NData+1)),ceil(sqrt(NData+1)),4);
imagesc(SumColorMap);
SumFinalrCellDensity=sum(FinalrCellDensity,2);
% caxis([min(SumFinalrCellDensity),max(SumFinalrCellDensity)])
dd=SumColorMap.*rangeMask;
dd(dd==0)=nan;
caxis([min(min(dd)),max(max(dd))])
axis equal
colorbar

cd(OutDir)
print(['Data#',num2str(iSurr),'.pdf'],'-dpdf')

FinalOut(iSurr).IndividualColorMap=IndividualColorMap;
FinalOut(iSurr).SumColorMap=SumColorMap;
FinalOut(iSurr).FinalPosUM=FinalPosUM;
FinalOut(iSurr).RotCorrect=RotCorrect;

end



%Circular power spectrumのデータを追加する。

tic
h0 = waitbar(0, 'データ解析中...'); % waitbar作成
for iSurr=1:NSurr+1 %数字注意      

%tangential距離フィルタ内の細胞密度を極座標の波形として算出する。
Center = [length(Bin1), length(Bin2)] * 0.5;
NAngle = 360;

% グリッドの座標
[X, Y] = meshgrid(1:length(Bin2), 1:length(Bin1));  % Y: row, X: column
dX = X - Center(2); % 列方向
dY = -(Y - Center(1)); % 行方向 (上下反転)

% 極座標に変換（1回だけ）
[theta, rho] = cart2pol(dX, dY); 
theta_deg = rad2deg(theta);
theta_deg(theta_deg < 0) = theta_deg(theta_deg < 0) + 360;

% 範囲マスク（1回だけ）
rangeMask = rho >= RangeUM(1)/BinWidth & rho < RangeUM(2)/BinWidth;

for iData = 1:NData
        rCellDensity = nan(NAngle, 1);
        rCellDensity2 = nan(NAngle, 1);
        
        for iAngle = 1:NAngle
            TempAng = iAngle - 1;
            AngRange = [TempAng - 5, TempAng + 5];

            % 角度範囲のマスク
            angleMask = theta_deg >= AngRange(1) & theta_deg < AngRange(2);
            fullMask = angleMask & rangeMask;

            % 適用と平均値算出
            Temp7 = FinalOut(iSurr).IndividualColorMap(:,:,iData) .* double(fullMask);
            Temp7(Temp7 == 0) = NaN;
            rCellDensity(iAngle) = nanmean(Temp7(:));
            
             Temp8 = FinalOut(iSurr).SumColorMap .* double(fullMask);
             Temp8(Temp8 == 0) = NaN;
             rCellDensity2(iAngle) = nanmean(Temp8(:));
            
        end

        FinalOut(iSurr).IndividualrCellDensity(:,iData) = rCellDensity;
        FinalOut(iSurr).SumrCellDensity = rCellDensity2;
        
        %Circular power spectrum analysis

        power_spec=[];
        theta_deg2 = 0:359;
        theta_rad = 2 * pi * theta_deg2 / 360;
        max_harmonic = 10;  % 1-10回の対称性まで調べる

        y = rCellDensity';
        N = length(y);

        for k = 1:max_harmonic
            A_n = sum(y .* exp(-1i * k * theta_rad)) / N;
            power_spec(k) = abs(A_n).^2;
        end

        FinalOut(iSurr).IndividualPowerSpect(:,iData)=power_spec;
        
        
        power_spec=[];
        theta_deg2 = 0:359;
        theta_rad = 2 * pi * theta_deg2 / 360;
        max_harmonic = 10;  % 1-10回の対称性まで調べる

        y = rCellDensity2';
        N = length(y);

        for k = 1:max_harmonic
            A_n = sum(y .* exp(-1i * k * theta_rad)) / N;
            power_spec(k) = abs(A_n).^2;
        end
        
        FinalOut(iSurr).SumPowerSpect(:,1)=power_spec;
        
end

d=FinalOut(iSurr).IndividualPowerSpect;
d1=d;
d2=d1./sum(d1);

FinalOut(iSurr).RatioIndividualPowerSpect=d2;
FinalOut(iSurr).AveRatioIndividualPowerSpect=sum(d2,2)/2;

waitbar(iSurr/(NSurr+1),h0);
toc

end
close(h0)

for i=1:NSurr+1
    RelatHexaPowerSpect(i)=FinalOut(i).AveRatioIndividualPowerSpect(6);
end

figure;
hold on;
for i=1:NSurr
    plot(2:2:10,FinalOut(i+1).AveRatioIndividualPowerSpect(2:2:10),'-','Color',0.8*[1 1 1]);
end
plot(2:2:10,FinalOut(1).AveRatioIndividualPowerSpect(2:2:10),'k-','LineWidth',3);
xlabel('Fold symmetry','FontSize',15)
ylabel('Power fraction','FontSize',15)


%% 高速化関数群

function result = applyRotationVectorized(posUM, rotMtx1, rotMtx2)
    % ベクトル化された座標回転
    if isempty(posUM)
        result = posUM;
        return;
    end
    
    % 第1回転 (XY平面)
    temp1 = rotMtx1 * posUM(:,1:2)';
    
    % 第2回転 (XZ平面)
    temp2 = rotMtx2 * [temp1(1,:); posUM(:,3)'];
    
    result = [temp2(1,:)', temp1(2,:)', temp2(2,:)'];
end

function unitVec3 = precomputeUnitVectors(angle1, angle2, unitVec)
    % 単位ベクトルの事前計算
    nAngles = length(angle1) * length(angle2);
    unitVec3 = zeros(nAngles, 3);
    
    m = 1;
    for i = 1:length(angle1)
        for j = 1:length(angle2)
            theta = deg2rad(angle1(i));
            phi = deg2rad(angle2(j));
            
            rMtx1 = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            rMtx2 = [cos(phi) -sin(phi); sin(phi) cos(phi)];
            
            unitVec2 = unitVec;
            temp22 = rMtx2 * [unitVec(1); unitVec(3)];
            unitVec2(1) = temp22(1);
            unitVec2(3) = temp22(2);
            temp11 = rMtx1 * [unitVec2(1); unitVec2(2)];
            unitVec2(1) = temp11(1);
            unitVec2(2) = temp11(2);
            
            unitVec3(m,:) = unitVec2;
            m = m + 1;
        end
    end
end

function exrPosUM = computeRelativePositionsVectorized(tempPosUM1, tempPosUM2, zlimMax)
    % 相対座標計算の高速化
    if isempty(tempPosUM1) || isempty(tempPosUM2)
        exrPosUM = [];
        return;
    end
    
    n1 = size(tempPosUM1, 1);
    n2 = size(tempPosUM2, 1);
    
    % より効率的な相対座標計算
    % kronを使用してtempPosUM1の各行をn2回繰り返し、
    % repmatを使用してtempPosUM2全体をn1回繰り返す
    repeleTempPosUM1 = kron(tempPosUM1, ones(n2, 1));
    repTempPosUM2 = repmat(tempPosUM2, n1, 1);
    
    trPosUM = repTempPosUM2 - repeleTempPosUM1;
    
    % Z制限の適用
    zLim = [0, zlimMax];
    validIdx = abs(trPosUM(:,3)) > zLim(1) & abs(trPosUM(:,3)) < zLim(2);
    exrPosUM = trPosUM(validIdx, :);
end

function kernel = precomputeGaussianKernel(sigma, binWidth)
    % ガウシアンカーネルの事前計算
    kernelSize = ceil(3 * sigma / binWidth);
    [X, Y] = meshgrid(-kernelSize:kernelSize, -kernelSize:kernelSize);
    kernel = exp(-(X.^2 + Y.^2) / (2 * (sigma/binWidth)^2));
    kernel = kernel / sum(kernel(:));
end

function gaussMap = GaussConvOptimized(posData, bin1, bin2, kernel)
    % 高速化されたガウシアンフィルタ
    if isempty(posData)
        gaussMap = zeros(length(bin1), length(bin2));
        return;
    end
    
    % ヒストグラム作成
    [N, ~, ~] = histcounts2(posData(:,1), posData(:,2), bin1, bin2);
    
    % 畳み込み
    gaussMap = conv2(N, kernel, 'same');
end

function [polarGrid, rangeMask] = precomputePolarGrid(bin1, bin2, rangeUM, binWidth)
    % 極座標グリッドの事前計算
    center = [length(bin1), length(bin2)] * 0.5;
    [X, Y] = meshgrid(1:length(Bin2), 1:length(Bin1));
    dX = X - center(2);
    dY = -(Y - center(1));
    
    [theta, rho] = cart2pol(dX, dY);
    theta_deg = rad2deg(theta);
    theta_deg(theta_deg < 0) = theta_deg(theta_deg < 0) + 360;
    
    rangeMask = rho >= rangeUM(1)/binWidth & rho < rangeUM(2)/binWidth;
    
    polarGrid.theta_deg = theta_deg;
    polarGrid.rho = rho;
end

function rCellDensity = computePolarDensityVectorized(colorMap, polarGrid, rangeMask, nAngle)
    % ベクトル化された極座標密度計算
    rCellDensity = nan(nAngle, 1);
    
    for iAngle = 1:nAngle
        tempAng = iAngle - 1;
        angRange = [tempAng - 5, tempAng + 5];
        
        angleMask = polarGrid.theta_deg >= angRange(1) & polarGrid.theta_deg < angRange(2);
        fullMask = angleMask & rangeMask;
        
        temp = colorMap .* double(fullMask);
        temp(temp == 0) = NaN;
        rCellDensity(iAngle) = nanmean(temp(:));
    end
end

function powerSpec = computeCircularPowerSpectrum(data, maxHarmonic)
    % 高速化されたフーリエ変換
    theta_deg = 0:359;
    theta_rad = 2 * pi * theta_deg / 360;
    N = length(data);
    
    powerSpec = zeros(1, maxHarmonic);
    
    for k = 1:maxHarmonic
        A_n = sum(data .* exp(-1i * k * theta_rad)) / N;
        powerSpec(k) = abs(A_n)^2;
    end
end