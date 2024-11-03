%% separate the calcium spikes into each trial
%load twoP2.mat     %twoP and twoP2 variables are loaded to workspace
deNoiseTraceF_F = twoP2.F_add.deltaF_F;   % column: loop number; row: n of neurons
[loopNum,cellNum] = size(deNoiseTraceF_F); %get the loop number and cell number by the denoised matrix
matTraceF_FwithTimestamp = zeros(loopNum,3,cellNum);
for i = 1:cellNum
    matTraceF_FwithTimestamp(:,3,i) = deNoiseTraceF_F(:,i);
end   %resize the matrix into row(loop),column(value+timestamp),z(neurons) pattern
%%
%load the timestamp data
disp('---xls.file, timestamp of current plane---');
[Temp1,Dir]=uigetfile('*.txt');
cd(Dir);
Time=importdata(Temp1);
FrameTime=(max(Time(:,1))-Time(1,1))/(length(Time)-1);
Time1=(Time(1,1):FrameTime:max(Time(:,1)))';  %adjust the frame time and make the step even
t=Time1;
t0=t-t(1);

%% load the behavioral .mat data of this training session
disp('---select the behavioral .mat file for the current session---');
[h,Dir]=uigetfile('*.mat');
cd(Dir);
load(h);

%%  get the time lag from imaging start to behavioral program start

%balance the difference between totalTime and sum of mLatency
timeLag = totalTime - sum(mLatency(:));
%change the matrix into the row-wise sum, and transpose into column vector
lat = mLatency';
lat = reshape(lat,[],1);
lat = cumsum(lat);
totalLatency = reshape(lat,[],160);  
totalLatency = totalLatency + timeLag/2;  
%totalLatency:   row1 -- end of cue
%                row2 -- end of post-cue 
%                row3 -- end of stim
%                row4 -- end of RW
%                row5 -- end of ITI
%                row6 -- time lapse to next trial


%% the latency between totalTime and sum of the mLatency is balanced

disp('------please input the latency between imaging starts and trial starts----------')

diffImaging = input("Imaging-behavioral start time lapse: ");   %usually imaging start earlier, 
                                                                % if later please input negative value


totalLatency = totalLatency + diffImaging; %adjust the time stamp of behavioral timestamp


%% separate the timestamp(t0) using the behavioral timestamp information 

vec1 = t0(:,1);                               %input of selected plane
thresholds = totalLatency(6,:);               %timepoint of every trial ended
thresholds1 = [totalLatency(1,1), thresholds,inf]; %the vector of the start point/each trial end/inf
trialLabels = discretize(vec1,thresholds1);   %important function discretize
t0(:,2) = trialLabels;
for i = 1:cellNum
 matTraceF_FwithTimestamp(:,1,i) = t0(:,2);  %first row, trial label
 matTraceF_FwithTimestamp(:,2,i) = t0(:,1);  %second row, timestamp
end
%so far the structure in matTraceF_FwithTimestamp is
%[trialNum,timeStamp,F_F] * z(neuronID)

%% separate the matrices for individual trials

matTrialLabels = unique(trialLabels);
trialMatrices = cell(length(matTrialLabels),1); %make the cells for individual trials
neuronMatrices = cell(1,cellNum);

for j = 1:cellNum
  for i = 1:trialNum    % 1-160 trials
     currentTrial = matTrialLabels(i);
     trialMatrices{i} = matTraceF_FwithTimestamp(matTraceF_FwithTimestamp(:,1) == currentTrial,:,j);
     neuronMatrices{i,j} = trialMatrices{i}; 
  end
end

%% separate the trials by contrast informat

neuronMatricesOnContrast = [num2cell(h.data1(:,1)),...   %trial#
                            num2cell(h.data1(:,2)),...   %trial result
                            num2cell(h.data1(:,7)),...   %trial contrast
                            neuronMatrices];             %neuronal activity info

neuronMatricesOnContrast = sortrows(neuronMatricesOnContrast,...
                                    [3 2],{'descend' 'ascend'});  
                           %sort the cell array by decending contrast and
                           %ascending result flag 

%% plot the neuronal activity in individual trials based on different contrast

%separate the 

waitPlot = waitbar(0,"plotting...");
%title("#1 to 10 neurons activity change in all 100% contrast trials")
for i = 1:cellNum  %first 9 neurons 
  if mod(i,10) == 1 
      figure;
  end
  waitbar(i/cellNum,waitPlot,"plotting...")
  subplotLabel = mod(i,10);
  if subplotLabel == 0
      subplotLabel = 10;
  end
  subplot(5,2,subplotLabel)
  title(append('ROI #',num2str(i)))
  cellLabel = i+3;
  for j = 1:40   %all 100% trials
      hold on
      %xplot(:,j) = contrastStruct.hundred{j,4}(:,2) - totalLatency(3,j);  
      %yplot(:,j) = contrastStruct.hundred{j,4}(:,3);
    plot(contrastStruct.hundred{j,cellLabel}(:,2) - totalLatency(3,j), ...
        contrastStruct.hundred{j,cellLabel}(:,3), ...
        '-','Color',[0 0 0 0.1],'LineWidth',2);
    xlim([-1 8])
    xticks(-1:1:8)
    xline(1,'--','Response Window')
    xline(5,'--')
    xlabel('s')
    ylabel('F_F')
    hold off
  end
end
close(waitPlot)

%%

%waitPlot = waitbar(0,"plotting...");
figure;
seismic = mymap('seismic');    %import colormap from mymap function ('help mymap' for more information)
colormap(seismic);
colorBar = colorbar;
colorBar.Label.String = 'Z score';
xlim([-1 8])
ylim([1 160])
xticks(-1:1:8)
xline(1,'--','Response Window')
xline(5,'--')
xlabel('s')
ylabel('z score')
waitPlot = waitbar(0,"plotting...");
set(gca,'YDir','reverse');
for i = 1:160   %trial number
    waitbar(i/160,waitPlot,"plotting...")
    hold on
    xplot = neuronMatricesOnContrast{i,4}(:,2) - totalLatency(3,i);
    yIncrement = 1;
    yPosition = 1 + i * yIncrement * ones(size(xplot));
    zScore = zscore(neuronMatricesOnContrast{i,4}(:,3));
    scatter(xplot,yPosition,30,zScore,'square','filled');  
end
hold off
close(waitPlot)


