%% summarize and analyze the data from 3 mice
disp('----------choose the corrAndDist.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M3 = corrAndDist;
%delete all the odd rows
CorrAndDist_M3(1:2:end,:) = [];
clearvars CorrAndDist
%
disp('----------choose the corrAndDist.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M9 = corrAndDist;
%delete all the odd rows
CorrAndDist_M9(1:2:end,:) = [];
clearvars CorrAndDist
%
disp('----------choose the corrAndDist.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('corrAndD*.mat');
cd(path)
load(fileName)
CorrAndDist_M17 = corrAndDist;
%delete all the odd rows
CorrAndDist_M17(1:2:end,:) = [];
clearvars CorrAndDist

% cat the data
corrAndDist_all = [CorrAndDist_M3; CorrAndDist_M9; CorrAndDist_M17];

%exclude the same ROI pair
corrAndDist_all = corrAndDist_all(corrAndDist_all(:,2) ~= 0,:);

%exclude the cell pair with more than 0.5 CC
corrAndDist_all = corrAndDist_all(corrAndDist_all(:,1) <= 0.5,:);

%sort the data based on the second column(untransformed XY distance)
corrAndDist_all = sortrows(corrAndDist_all,2);

%exclude the cell pair with xyz distance less than 30
corrAndDist_all1 = corrAndDist_all(corrAndDist_all(:,3) > 30,:);

% 
figure;
hold on
xlim([0 50]);
ylim([-0.1 0.2]);
scatter(corrAndDist_all1(:,2),corrAndDist_all1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('CC and untransformed Distance of ROIs in 3 mice'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndDist_all1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndDist_all1(:,2) >= distanceBin(i) & corrAndDist_all1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndDist_all1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(2:end)/2,corrBin,'k','LineWidth',2);  %plot the average correlation coefficient of the distance bin
xlabel('Distance between the paired ROI centroids');     %add the x-axis label
ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
set(gca,'FontSize',16)

%% summarize and analyze the data from 3 mice
disp('----------choose the corrAndTransDist.mat of the 1st animal(M3)----------')
[fileName,path] = uigetfile('corrAndTrans*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M3 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M3(1:2:end,:) = [];
clearvars CorrAndTransDist
%
disp('----------choose the corrAndTransDist.mat of the 2nd animal(M9)----------')
[fileName,path] = uigetfile('corrAndTransDist*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M9 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M9(1:2:end,:) = [];
clearvars CorrAndTransDist
%
disp('----------choose the corrAndTransDist.mat of the 3rd animal(M17)----------')
[fileName,path] = uigetfile('corrAndTransDist*.mat');
cd(path)
load(fileName)
CorrAndTransDist_M17 = corrAndTransDist;
%delete all the odd rows
CorrAndTransDist_M17(1:2:end,:) = [];
clearvars CorrAndTransDist

% cat the data
corrAndTransDist_all = [CorrAndTransDist_M3; CorrAndTransDist_M9; CorrAndTransDist_M17];

%exclude the same ROI pair
corrAndTransDist_all = corrAndTransDist_all(corrAndTransDist_all(:,2) ~= 0,:);

%exclude the cell pair with 1 correlation coefficient
corrAndTransDist_all = corrAndTransDist_all(corrAndTransDist_all(:,1) <= 0.5,:);
%sort the data based on the second column(transDist)
corrAndTransDist_all = sortrows(corrAndTransDist_all,2);

%exclude the cell pair with xyz distance less than 30
corrAndTransDist_all1 = corrAndTransDist_all(corrAndTransDist_all(:,3) > 30,:);

figure;
hold on
xlim([0 50]);
ylim([-0.05 0.1]);

%scatter(corrAndTransDist_all1(:,2),corrAndTransDist_all1(:,1),10,[0.8 0.8 0.8],'Filled'); %plot the scatter plot of correlation coefficient and distance
title('CC and tangential Distance of ROIs in 3 mice'); %add the title
set(gca,'FontSize',16)
%plot the average correlation coefficient of ROIs with the different distance bin
distanceBin = 0:5:max(corrAndTransDist_all1(:,2));       %create the distance bin
corrBin = zeros(length(distanceBin)-1,1);       %create a vector to store the average correlation coefficient of each distance bin
for i = 1 : length(distanceBin)-1
    tempIdx = find(corrAndTransDist_all1(:,2) >= distanceBin(i) & corrAndTransDist_all1(:,2) < distanceBin(i+1)); %find the index of the correlation coefficient in the distance bin
    corrBin(i) = mean(corrAndTransDist_all1(tempIdx,1));  %calculate the average correlation coefficient of the distance bin
end
plot(distanceBin(2:end)/2,corrBin,'Color',[0.1 0.7 0.1],'LineWidth',4);  %plot the average correlation coefficient of the distance bin
xlabel('Distance between the paired ROI centroids');     %add the x-axis label
ylabel('Correlation Coefficient of ROI pairs');          %add the y-axis label
set(gca,'FontSize',16)

