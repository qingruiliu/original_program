%% select the result data at same contrast leve;
%hit
disp('----------choose the Hit----------')
[fileName,path] = uigetfile('*Hit*.mat');
cd(path)
load(fileName)
Hit_temp = tempVar;

%Miss
disp('----------choose the Miss----------')
[fileName,path] = uigetfile('*Miss*.mat');
cd(path)
load(fileName)
Miss_temp = tempVar;

%FA
disp('----------choose the FA----------')
[fileName,path] = uigetfile('*FA*.mat');
cd(path)
load(fileName)
FA_temp = tempVar;

%CR
disp('----------choose the CR----------')
[fileName,path] = uigetfile('*CR*.mat');
cd(path)
load(fileName)
CR_temp = tempVar;

% plot the mean of real data and the mean of 1000 surrogates for different result on the same figure
% Hit
sizeHit = length(Hit_temp);
sizeMiss = length(Miss_temp);
sizeFA = length(FA_temp);
sizeCR = length(CR_temp);

allCorr_In20 = [Hit_temp(:,1);Miss_temp(:,1);FA_temp(:,1);CR_temp(:,1)];

surrogateTime = 1000;
hitSurroData = zeros(sizeHit,surrogateTime);
missSurroData = zeros(sizeMiss,surrogateTime);
faSurroData = zeros(sizeFA,surrogateTime);
crSurroData = zeros(sizeCR,surrogateTime);

%make the random surrogate data for 1000 times and make the grouped scatter plot
hitColor = [0.25 0.8 0.25];
missColor = [1 0.54 0.1];
faColor = [0.83 0.14 0.14];
crColor = [0.27 0.25 0.8];
wb = waitbar(0,'Please wait...');
figure
xlim([0 5])
xticks([1 2 3 4])
xticklabels({'Hit','Miss','FA','CR'})

hold on
for i = 1:surrogateTime
    waitbar(i/surrogateTime,wb);
    %make the surrogate hit data from allCorr_In20
    randHitIdx = randi([1,length(allCorr_In20)],sizeHit,1);
    randMissIdx = randi([1,length(allCorr_In20)],sizeMiss,1);
    randFAIdx = randi([1,length(allCorr_In20)],sizeFA,1);
    randCRIdx = randi([1,length(allCorr_In20)],sizeCR,1);

    hitSurroData(:,i) = allCorr_In20(randHitIdx);
    missSurroData(:,i) = allCorr_In20(randMissIdx);
    faSurroData(:,i) = allCorr_In20(randFAIdx);
    crSurroData(:,i) = allCorr_In20(randCRIdx);

    %plot the mean of surrogate data in each loop
   
end
close(wb)

%plot the mean of each surrogate data
scatter(ones(sizeHit,1),mean(hitSurroData,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(2*ones(sizeMiss,1),mean(missSurroData,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(3*ones(sizeFA,1),mean(faSurroData,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(4*ones(sizeCR,1),mean(crSurroData,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)

%plot the mean of real data as a bar
plot(1,mean(Hit_temp(:,1)),'o','MarkerSize',8,'MarkerFaceColor',hitColor,'MarkerEdgeColor',hitColor)
plot(2,mean(Miss_temp(:,1)),'o','MarkerSize',8,'MarkerFaceColor',missColor,'MarkerEdgeColor',missColor)
plot(3,mean(FA_temp(:,1)),'o','MarkerSize',8,'MarkerFaceColor',faColor,'MarkerEdgeColor',faColor)
plot(4,mean(CR_temp(:,1)),'o','MarkerSize',8,'MarkerFaceColor',crColor,'MarkerEdgeColor',crColor)

ylim([0.1 0.2])
yticks([0.1:0.05:0.2])
set(gca,'FontSize',16)
%set the tick length as 0
set(gca,'TickLength',[0 0])
ylabel('Average correlation')

%% save the figure as .eps file
print('01_resultCompare', '-depsc', '-painters');



%% select the result data at same contrast leve;
%hit
disp('----------choose the 100----------')
[fileName,path] = uigetfile('*_100_*.mat');
cd(path)
load(fileName)
temp_100 = tempVar;

%Miss
disp('----------choose the 10----------')
[fileName,path] = uigetfile('*_10_*.mat');
cd(path)
load(fileName)
temp_10 = tempVar;

%FA
disp('----------choose the 1----------')
[fileName,path] = uigetfile('*_1_*.mat');
cd(path)
load(fileName)
temp_1 = tempVar;

%CR
disp('----------choose the 0.1----------')
[fileName,path] = uigetfile('*_01_*.mat');
cd(path)
load(fileName)
temp_01 = tempVar;

%% plot the mean of real data and the mean of 1000 surrogates for different result on the same figure
% Hit
size100 = length(temp_100);
size10 = length(temp_10);
size1 = length(temp_1);
size01 = length(temp_01);

allCorr_In10 = [temp_100(:,1);temp_10(:,1);temp_1(:,1);temp_01(:,1)];

surrogateTime = 1000;
surroData100 = zeros(size100,surrogateTime);
surroData10 = zeros(size10,surrogateTime);
surroData1 = zeros(size1,surrogateTime);
surroData01 = zeros(size01,surrogateTime);

%make the random surrogate data for 1000 times and make the grouped scatter plot
hitColor = [0.25 0.8 0.25];
missColor = [1 0.54 0.1];
faColor = [0.83 0.14 0.14];
crColor = [0.27 0.25 0.8];
wb = waitbar(0,'Please wait...');
figure
xlim([0 5])
xticks([1 2 3 4])
xticklabels({'100%','10%','1%','0.1%'})

hold on
for i = 1:surrogateTime
    waitbar(i/surrogateTime,wb);
    %make the surrogate hit data from allCorr_In20
    rand100Idx = randi([1,length(allCorr_In10)],size100,1);
    rand10Idx = randi([1,length(allCorr_In10)],size10,1);
    rand1Idx = randi([1,length(allCorr_In10)],size1,1);
    rand01Idx = randi([1,length(allCorr_In10)],size01,1);

    surroData100(:,i) = allCorr_In10(rand100Idx);
    surroData10(:,i) = allCorr_In10(rand10Idx);
    surroData1(:,i) = allCorr_In10(rand1Idx);
    surroData01(:,i) = allCorr_In10(rand01Idx);

    %plot the mean of surrogate data in each loop
   
end
close(wb)

%plot the mean of each surrogate data
scatter(ones(size100,1),mean(surroData100,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(2*ones(size10,1),mean(surroData10,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(3*ones(size1,1),mean(surroData1,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)
scatter(4*ones(size01,1),mean(surroData01,2),10,[0.7 0.7 0.7],'filled','jitter','on','jitterAmount',0.1)

%plot the mean of real data as a bar
scatter(1,mean(temp_100(:,1)),80,crColor,'filled')
scatter(2,mean(temp_10(:,1)),80,crColor,'filled') %'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8
scatter(3,mean(temp_1(:,1)),80,crColor,'filled')
scatter(4,mean(temp_01(:,1)),80,crColor,'filled')
ylim([0.1 0.22])
yticks(0.1:0.04:0.24)
set(gca,'FontSize',16)
%set the tick length as 0
set(gca,'TickLength',[0 0])
ylabel('Average correlation')

%% save the figure as .eps file
print('low_resultCompare', '-depsc', '-painters');


