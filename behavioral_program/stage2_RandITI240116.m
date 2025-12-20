%% modified program from stageTwoRasterPlot.m 24.1.5
% new functions (original)
%1. randomized post-cue period following normal distribution 
% with mean=1, SD = 0.1
%2. modified ploting tic/toc timer setting
% 
% Latest improvements:
%1. 优化视觉刺激流畅度 - 参考visual_stimulation_gui的实现方式
%2. 移除in/out response window统计功能 - 简化GUI和代码
%3. 增加lick rate monitor - 在ITI最后1秒监测小鼠lick rate并保存数据
timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor, h the gray color background
global h 
PsychDefaultSetup(2);

% 添加性能优化的Screen偏好设置
Screen('Preference', 'ConserveVRAM', 4096); % 优化显存使用
Screen('Preference', 'VBLTimestampingMode', 4); % 高精度VBL时间戳
Screen('Preference', 'SkipSyncTests', 0); % 确保同步测试通过
Screen('Preference', 'VisualDebugLevel', 0); % 关闭视觉调试信息
Screen('Preference', 'SuppressAllWarnings', 1); % 抑制警告信息

Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);
h.screenNumber = max(Screen('Screens'));
h.white = WhiteIndex(h.screenNumber);
h.grey = h.white / 2;

% 优化窗口打开参数
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution); 
h.ifi = Screen('GetFlipInterval',h.window); 

% 设置最高优先级和实时调度
h.topPriorityLevel = MaxPriority(h.window); 
Priority(h.topPriorityLevel);

% 预热显卡和优化缓存
for i = 1:5
    Screen('Flip', h.window);
end

%% initialize sound configuration
InitializePsychSound;

%open psych-audio port
h.sampleF = 48000;
% Find the device ID for the audio interface
deviceList = PsychPortAudio('GetDevices');

% find the device ID for the Steinberg UR12
device_id = [];
for i = 1:length(deviceList)
    if startsWith(deviceList(i).DeviceName, 'Steinberg UR12')
        device_id = i;
        disp(['Steinberg UR12 found with Device ID: ' num2str(device_id)]);
        break;
    end
end
if isempty(device_id)
    error('Steinberg UR12 device not found.');
end

h.audioHandle = PsychPortAudio('Open',device_id-1, 1, 1, h.sampleF, 2);  %use PsychPortAudio('GetDevices') to find steinberg UR12, and change the first number with UR12 number
PsychPortAudio('Volume', h.audioHandle, 0.1);      %auditory cue volume

%pre-allocate audio buffer
[myBeep, samplingRate] = MakeBeep(10000, 0.1, h.sampleF);
buffer = [myBeep;myBeep];
PsychPortAudio('FillBuffer',h.audioHandle,buffer);

%% start communication
h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
h.sensorPin = 'D13';
h.waterPumpPin = 'D9';

%% Gabor presetting - 优化设置，参考visual_stimulation_gui
%使用与visual_stimulation_gui相同的设置方式
[h.width, h.height] = Screen('WindowSize', h.window);
h.gaborDimPix = max(h.width, h.height); % 使用屏幕尺寸

%other parameters - 与visual_stimulation_gui保持一致
h.contrast = 1.0;
h.phase = 0;
h.orientationTarget = 0; %vertical
h.orientationNontarget = 90; %horizontal

%使用与visual_stimulation_gui相同的空间频率设置
h.numCycles = 7; % 与visual_stimulation_gui相同
h.freq = h.numCycles / h.gaborDimPix; % 与visual_stimulation_gui相同的计算方式

%使用与visual_stimulation_gui相同的参数
h.sigma = h.gaborDimPix;
h.aspectRatio = 1;
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;

%创建Gabor纹理，与visual_stimulation_gui保持一致
h.gratingtex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

%属性矩阵与visual_stimulation_gui相同
h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];

% 优化的相位增量和帧设置
h.phasePerFrame = 5 * pi; % 与visual_stimulation_gui完全相同
h.waitframes = 1; % 与visual_stimulation_gui相同

% 预计算一些常用值以提高性能
h.halfIfi = 0.5 * h.ifi;
h.frameAdvance = (h.waitframes - 0.5) * h.ifi;

%% program variables
trialNumLimit = 300;
h.licktrial = 1;
lickTrialLimit = 150;
h.lickdata = {};
h.data1 = zeros(trialNumLimit,4);
h.totalLickTimes = 0;

% Lick count monitoring variables - 热图版本（高时间分辨率）
% 时间窗选择说明：
% - 50ms: 极高分辨率，但可能噪声较大，适合精细分析
% - 100ms: 高分辨率与稳定性的平衡，推荐设置
% - 200ms+: 更平滑但分辨率较低
% 注意: 显示的是lick count (times/trial/100ms)，而非lick rate (Hz)
h.lickRateTimeWindow = 0.1; % 时间窗大小（秒）- 100ms窗口
h.lickRateTimeStep = 0.05; % 时间步长（秒）- 50ms步长，提供更平滑的热图
h.trialLickTimes = []; % 当前trial中所有lick的时间戳
h.lickRateTimeCourse = {}; % 保存每个trial的lick count时间过程
h.timeAxis = 0:h.lickRateTimeStep:11; % 时间轴（0-11秒）

% 热图数据矩阵 - 行为trial，列为时间点
h.lickRateMatrix = zeros(trialNumLimit, length(h.timeAxis));
h.heatmapHandle = []; % 热图句柄

%create random ITI time length ranging from 4 ~ 6s.
h.ITIperiod = 4 + rand([1 trialNumLimit])*2;
%% display UI, waiting for the initialization
infoUI();

%% counterUI - 移除in/out RW统计，简化GUI
screenSize = get(0,'Screensize'); screenSize(3) = screenSize(3)/2;       
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);          %open the lick monitor UI

%total counter UI - 调整为一行显示，减少高度
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.05 0.95 0.08 0.03],'String','Total Trial','BackgroundColor',[0 1 1],...
    'FontSize',12);
 h.totalTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.05 0.92 0.08 0.025],'FontSize',12); 
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.15 0.95 0.08 0.03],'String','Lick Trial','BackgroundColor',[0 1 1],...
    'FontSize',12);
 h.lickTrialNumUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.15 0.92 0.08 0.025],'FontSize',12); 
 % function of hit trial number
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.25 0.95 0.08 0.03],'String','Hit Trial','BackgroundColor',[1 1 0],...
    'FontSize',12);
 h.hitTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.25 0.92 0.08 0.025],'FontSize',12); 

 %visual stimulation parameters - 调整位置
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.35 0.95 0.06 0.03],'String','VS','BackgroundColor',[1 1 1],...
     'FontSize',12);
 h.tempF =uicontrol(f,'Style','edit','String','1','Units','normalized',...
    'Position',[0.35 0.92 0.06 0.025],'FontSize',12); 
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.43 0.95 0.06 0.03],'String','RW','BackgroundColor',[1 1 1],...
     'FontSize',12);
 h.spatF =uicontrol(f,'Style','edit','String','4','Units','normalized',...
    'Position',[0.43 0.92 0.06 0.025],'FontSize',12); 

 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.51 0.95 0.06 0.03],'String','ITI','BackgroundColor',[1 1 1],...
     'FontSize',12);
 h.duration =uicontrol(f,'Style','edit','String','4','Units','normalized',...
    'Position',[0.51 0.92 0.06 0.025],'FontSize',12); 

% trial raster - 调整位置以适应新的GUI布局
h.trialRaster = axes(f,'Position',[0.1 0.5 0.3 0.35],'FontSize',14);
title('trial raster');
xlabel('seconds');
ylabel('trial number');
xlim([0 11]);
ylim([0 300]);
set(h.trialRaster,'Ydir','reverse');
hold on

% lick rate monitor - 热图显示版本
h.lickRateMonitor = axes(f,'Position',[0.55 0.5 0.3 0.35],'FontSize',14);
title('Lick Count Heatmap (times/trial/100ms)');
xlabel('Time in trial (s)');     
xlim([0 11]); % 与trial raster相同的时间范围
ylabel('Trial number');
ylim([0 300]); % 初始Y轴范围设为300，与trial raster一致
set(h.lickRateMonitor,'YDir','reverse'); % 从上到下排列，与trial raster一致
hold on

%open the back camera
% Get camera resolutions
h.backCam = webcam(1);
h.frontCam = webcam(2);

backRes = str2double(strsplit(h.backCam.Resolution,'x'));
frontRes = str2double(strsplit(h.frontCam.Resolution,'x'));

% Calculate aspect ratios
backAspect = backRes(1) / backRes(2);
frontAspect = frontRes(1) / frontRes(2);

% Set UI positions based on aspect ratios (normalized units) - 调整位置以适应新布局
backHeight = 0.28;
backWidth = backAspect * backHeight;
frontHeight = 0.28;
frontWidth = frontAspect * frontHeight;

% Place back camera UI - 调整位置
h.backCamUI = axes(f, 'Position', [0.45 0.1 backWidth backHeight]);
h.im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', h.backCamUI);
preview(h.backCam, h.im);
text(h.backCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 12);  % Add label below the image

% Place front camera UI - 调整位置
h.frontCamUI = axes(f, 'Position', [0.05 0.1 frontWidth frontHeight]);
title(h.frontCamUI, 'Back Camera');
h.im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', h.frontCamUI);
preview(h.frontCam, h.im2);
text(h.frontCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 12);  % Add label below the image

%% define the timers
h.tLickCounter = timer('ExecutionMode', 'fixedRate', 'Period', 0.01,...
                             'TimerFcn',@(src,event)pinStatusChanged);
%global lick report timer with 10ms period
h.tRefractory = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.1,...
    'StartFcn',@(src,event)refStart,...
    'TimerFcn',@(src,event)refEnd);
%% 3 minutes countdown 
countDown;

%% start tic 
totalTic = tic;
h.inOrOutRW = [];
h.inRWCounter = 0;
h.lickInRWOneTrial = 0;
h.hitTrialNumber = 0;
h.earlyTrialNumber = 0;
h.trialNum = [];
h.rwTime = [];
h.trialLickTimes = []; % 初始化trial内lick时间记录
mLatency = zeros(trialNumLimit,6);   %latency test
%% main loop
start(h.tLickCounter);   %start the lick counter
for trialNum = 1:1000
       h.trialNum = trialNum;
    % 安全停止refractory timer
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end                                             %stop the tRefractory
    
    % 确保tLickCounter在trial间隙是停止的，避免误触发
    if strcmp(h.tLickCounter.Running,'on')
        stop(h.tLickCounter);
    end
    
h.lickInRWOneTrial = 0;
h.trialLickTimes = []; % 重置当前trial的lick时间记录

    if trialNum <= trialNumLimit  % within 300 trials (about 50 mins)
    h.trialGlobalTic = tic;
    disp('------------------new trial-----------------------')
    fprintf('trialNum = %s\n',num2str(trialNum))
    set(h.totalTrialNumUI,'String',num2str(trialNum));
   
     if trialNum > 1 
     mLatency(trialNum,6) = toc(latencyTic);
     end

    latencyTic = tic;
    trialCue;      %play the auditory cue of the single trial
    mLatency(trialNum,1) = toc(latencyTic);

    h.postCuePeriodlatencyTic = tic;
    postCuePeriod; %1 second after the auditory cue is given
    mLatency(trialNum,2) = toc(h.postCuePeriodlatencyTic);

    latencyTic = tic;
    visiStim;      %present the visual stimulation
    mLatency(trialNum,3) = toc(latencyTic);

    latencyTic = tic;
    responseWindow;
    mLatency(trialNum,4) = toc(latencyTic);

    latencyTic =tic;
    ITIperiod;      %intertrial interval
    mLatency(trialNum,5) = toc(latencyTic);
    
    latencyTic =tic;

    % 计算并绘制当前trial的lick rate时间过程
    calculateAndPlotLickRate(trialNum);

    fprintf('In this trial, lick in RW = %s times! \n',num2str(h.lickInRWOneTrial));
    h.data1(trialNum,1) = trialNum;
    h.data1(trialNum,2) = h.lickInRWOneTrial;
    h.data1(trialNum,3) = length(h.trialLickTimes); % 保存总lick次数
    h.data1(trialNum,4) = toc(h.trialGlobalTic); %time length for individual trials
    else
        break
    end
end
stop(h.tLickCounter);
delete(h.tLickCounter);
delete(h.tRefractory);
disp('---------------------finished: 50 licked trials!---------------')
totalTime = toc(totalTic);
totalLickTimes = 0;
for i = 1 : numel(h.lickdata)
    totalLickTimes = totalLickTimes + length(h.lickdata{i});
end
PsychPortAudio('Close',h.audioHandle); 
sca;
fprintf('>> total time cost:  %s minutes %s seconds \n',num2str(floor((totalTime)/60)),num2str(mod(totalTime,60)));
fprintf('>> total lick times: %s times \n',num2str(h.inRWCounter));
h.lickdata;
h.data1;
h.lickRateTimeCourse; % 显示lick rate时间过程数据
h.lickRateMatrix; % 显示lick rate热图矩阵数据
%将数据保存到文件，命名为日期+mouseID+conditioning日期
save([datestr(now,'yyyymmdd'),'_',h.mouseID{1},'_cond.mat'], 'h');

%% functions
function calculateAndPlotLickRate(trialNum)
global h
    % 计算滑动窗口内的lick count (times/100ms)
    lickCount = zeros(size(h.timeAxis));
    
    for i = 1:length(h.timeAxis)
        windowStart = h.timeAxis(i) - h.lickRateTimeWindow/2;
        windowEnd = h.timeAxis(i) + h.lickRateTimeWindow/2;
        
        % 计算在当前时间窗内的lick次数（直接计数，不除以时间窗）
        licksInWindow = sum(h.trialLickTimes >= windowStart & h.trialLickTimes <= windowEnd);
        lickCount(i) = licksInWindow; % lick times per 100ms window
    end
    
    % 保存当前trial的lick count时间过程
    h.lickRateTimeCourse{trialNum} = lickCount;
    h.lickRateMatrix(trialNum, :) = lickCount;
    
    % 使用热图显示
    updateLickRateHeatmap(trialNum);
end

function updateLickRateHeatmap(trialNum)
global h
    % 只显示已完成的trials
    dataToShow = h.lickRateMatrix(1:trialNum, :);
    
    % 如果是第一个trial，初始化热图
    if trialNum == 1
        h.heatmapHandle = imagesc(h.lickRateMonitor, h.timeAxis, 1:trialNum, dataToShow);
        colormap(h.lickRateMonitor, 'hot'); % 使用热色图
        colorbar(h.lickRateMonitor);
        
        % 设置坐标轴 - 从上到下排列，与trial raster一致
        set(h.lickRateMonitor, 'YDir', 'reverse'); % Y轴反向，从上到下
        xlabel(h.lickRateMonitor, 'Time in trial (s)');
        ylabel(h.lickRateMonitor, 'Trial number');
        title(h.lickRateMonitor, 'Lick Count Heatmap (times/trial/100ms)');
        
        % 设置轴范围 - 初始设为300
        xlim(h.lickRateMonitor, [0 11]);
        ylim(h.lickRateMonitor, [0.5 300.5]); % 与trial raster保持一致的300范围
        
        % 添加阶段标记线
        hold(h.lickRateMonitor, 'on');
        plot(h.lickRateMonitor, [1 1], [0.5 300.5], 'w--', 'LineWidth', 1); % Post-cue结束
        plot(h.lickRateMonitor, [2 2], [0.5 300.5], 'w--', 'LineWidth', 1); % VS结束
        plot(h.lickRateMonitor, [6 6], [0.5 300.5], 'w--', 'LineWidth', 1); % RW结束
        text(h.lickRateMonitor, 0.5, 10, 'Cue', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 1.5, 10, 'VS', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 4, 10, 'RW', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 8.5, 10, 'ITI', 'Color', 'white', 'FontSize', 8);
        hold(h.lickRateMonitor, 'off');
    else
        % 更新现有热图
        set(h.heatmapHandle, 'CData', dataToShow, 'YData', 1:trialNum);
        % Y轴范围固定为300，不再动态调整
        % ylim(h.lickRateMonitor, [0.5 300.5]); % 已在初始化时设置
    end
    
    % 设置colorbar范围 - 针对lick count调整
    caxis(h.lickRateMonitor, [0 3]); % 0-3 licks/100ms范围，适应lick count显示
    
    drawnow;
end

function adjustTimeWindow(newWindow, newStep)
% 动态调整时间窗参数的工具函数
% 用法: adjustTimeWindow(0.05, 0.025)  % 50ms窗口，25ms步长
%      adjustTimeWindow(0.1, 0.05)   % 100ms窗口，50ms步长
% 注意: 热图显示的是lick count (times/trial/window)，而非lick rate
global h
    if nargin >= 1
        h.lickRateTimeWindow = newWindow;
        fprintf('Time window adjusted to %.1f ms\n', newWindow*1000);
        fprintf('Heatmap now shows lick count per %.1f ms window\n', newWindow*1000);
    end
    if nargin >= 2
        h.lickRateTimeStep = newStep;
        h.timeAxis = 0:h.lickRateTimeStep:11;
        % 重新初始化矩阵尺寸
        h.lickRateMatrix = zeros(size(h.lickRateMatrix, 1), length(h.timeAxis));
        fprintf('Time step adjusted to %.1f ms\n', newStep*1000);
        fprintf('Matrix resized to %dx%d\n', size(h.lickRateMatrix, 1), size(h.lickRateMatrix, 2));
    end
end

function checkTimerStatus()
% 调试用函数：检查timer状态
global h
    fprintf('=== Timer Status Check ===\n');
    fprintf('tLickCounter: %s\n', h.tLickCounter.Running);
    fprintf('tRefractory: %s\n', h.tRefractory.Running);
    fprintf('Current phase: %d (postCue:-1, normal:0, RW:1, ITI:2)\n', h.inOrOutRW);
    fprintf('========================\n');
end

function infoUI(~,~)
 global h
 prompt = {'mouseID', 'trainStage', 'dayNumber', 'saveDir'};
 dlgtitle = 'mouse information';
 dims = [1 35];
 definput = {'', '', '', '/Users/liuqr/files/MATLAB相关/code ref/test code'};
 h.mouseID = inputdlg(prompt, dlgtitle, dims, definput);
end

function updateVbl(~,~)
 global h
 h.vbl = Screen('Flip',h.window);
 h.vblt0 = h.vbl;
end

function countDown(~,~)
 countdownDuration = 60;
 disp('---------Countdown started...check the mouse and lick spout!!!--------------------')

 for remainingSeconds = countdownDuration :-1 :0
     fprintf('Time remaining: %d seconds \n',remainingSeconds);
     pause(1);
 end
 disp('Countdown finished! Start conditioning!')
end

%% trial procedure functions
function trialCue(~,~)
global h
PsychPortAudio('Start',h.audioHandle,1,0,1);
[startTime, endPositionSecs, xruns, estStopTime] =PsychPortAudio('Stop',h.audioHandle,1,1);
disp('>>Trial Cue ended!! Post-cue period starts!!')
end

function postCuePeriod(~,~)
global h
    h.inOrOutRW = -1;  %inOrOutRW value: postCuePeriod -1, RW 1, Visi and ITI 0.
    
    % 安全启动lick counter，避免重复启动
    if ~strcmp(h.tLickCounter.Running,'on')
        start(h.tLickCounter); 
    end
    
    h.postCueTime = tic;
    while toc(h.postCueTime) <=  1     %fixed 1 seconds
        % 添加短暂等待，减少CPU占用，避免忙等
        pause(0.001);
    end
    disp('>>post-cue period finished! Visual stimulation starts!')
end

function visiStim(~,~)
  global h
  h.inOrOutRW = 0; 
  
  % 优化：清空绘图缓冲区并设置背景
  Screen('FillRect', h.window, h.grey);
  
  % 使用高精度VBL同步
  vbl = Screen('Flip', h.window);
  vblt0 = vbl;
  startTime = vbl;
  
  % 重置相位到初始值
  h.propertiesMat(1) = 0;
  
  VSLength = 1; % 视觉刺激持续1秒
  
  % 优化的刺激呈现循环
  frameCount = 0;
  while (vbl - vblt0) <= VSLength
    frameCount = frameCount + 1;
    
    % 高效的纹理绘制
    Screen('DrawTexture', h.window, h.gratingtex, [], [], h.orientationTarget, [], [], [], [],...
          kPsychDontDoRotation, h.propertiesMat');
    
    % 计算下一个精确的翻转时间
    nextFlipTime = startTime + frameCount * h.ifi;
    
    % 高精度翻转，使用预测的下一帧时间
    vbl = Screen('Flip', h.window, nextFlipTime - 0.5 * h.ifi);
    
    % 更新相位
    h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
  end
  
  % 刺激结束，清除屏幕并翻转
  Screen('FillRect', h.window, h.grey);
  Screen('Flip', h.window);
  
  h.inOrOutRW = 1;
  disp('>>Visual stimulation ended!! Response window starts!!')
end

function responseWindow(~,~)
global h
 h.rwTime = tic;
 %h.inOrOutRW = 1;  %set as true when RW starts  
  while toc(h.rwTime) <= 4    %the RW last for 4 seconds
  end
  if h.lickInRWOneTrial > 0
      h.licktrial = h.licktrial + 1;
  end
h.inOrOutRW = 0;   %reset the inOrOutRW
disp('>>Response window ended!! ITI start!')
end

function ITIperiod(~,~)
global h
    h.inOrOutRW = 2;
    h.ITITime = tic;
    
    % ITI等待
    itiDuration = h.ITIperiod(h.trialNum);
    while toc(h.ITITime) <= itiDuration
        WaitSecs(0.01); % 短暂等待，避免占用过多CPU
    end
    
    fprintf('>> ITI period finished!!! Time length is %s\n',...
        num2str(h.ITIperiod(h.trialNum)))
end

function pinStatusChanged(~,~)
global h
    % 安全检查：确保timer处于正确状态
    if ~strcmp(h.tLickCounter.Running,'on')
        % 如果lick counter没有运行，说明可能处于异常状态，直接返回
        return;
    end
    
    RWflag = h.inOrOutRW;
    lickFlag = true;
    trialFlag = true;
    lickTimesReporter = 0;
    pinValue = readDigitalPin(h.a,h.sensorPin);
    
    % 记录lick时间戳（相对于trial开始时间）
    if pinValue == true && lickFlag == true
        trialElapsedTime = toc(h.trialGlobalTic);
        h.trialLickTimes = [h.trialLickTimes, trialElapsedTime];
        
        switch RWflag
            case 1 % lick in RW
                h.inRWCounter = h.inRWCounter + 1;
                if h.lickInRWOneTrial < 1 %只有第一次lick被奖励
                    if round(toc(h.rwTime),3) > 10
                        plot(h.trialRaster,1.001,h.trialNum,'.g');   %avoid the larger than 10 problem after rounding
                    else
                        plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.g');  %the first lick in RW is marked as green dot
                    end
                   writeDigitalPin(h.a,'D9',1);
                   h.hitTrialNumber = h.hitTrialNumber + 1;
                   set(h.hitTrialNumUI,'String',num2str(h.hitTrialNumber));
                else
                    plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.','Color',[0.6 0.6 0.6]); %other licks in RW are marked as blue dot
                end
                fprintf('lick in RW %s ',num2str(h.inRWCounter));
                start(h.tRefractory);
                lickFlag = false;     
           if lickFlag == false
             if trialFlag == true
             trialFlag = false;
             fprintf('licktrial = %s ',num2str(h.licktrial));
             set(h.lickTrialNumUI,'String',num2str(h.licktrial));
             end
             lickTimesReporter = lickTimesReporter + 1;
             h.lickInRWOneTrial = h.lickInRWOneTrial + lickTimesReporter;
             tocReporter = toc(h.rwTime);
             fprintf('licktime = %s  @%s seconds \n',num2str(h.lickInRWOneTrial),num2str(tocReporter));
             h.lickdata{h.licktrial}(h.lickInRWOneTrial,1) = h.lickInRWOneTrial;
             h.lickdata{h.licktrial}(h.lickInRWOneTrial,2) = tocReporter;
           end
            case -1 % early lick reset the post-cue period timer
                    h.postCueTime = tic;
                    disp('!!Early lick! Reset post cue period timer!!')
                    start(h.tRefractory);
            case 2 % ITI期间的lick
                   plot(h.trialRaster,round(toc(h.ITITime),3)+5,h.trialNum,'.','Color',[0.6 0.6 0.6]); %plot the ITI licking in 5~9 s ITI window
                   start(h.tRefractory);
            otherwise % 其他时期的lick（如视觉刺激期间）
                   if RWflag == 0 && exist('h.vbl','var') && exist('h.vblt0','var')
                       % 在视觉刺激窗口期间的早期lick
                       plot(h.trialRaster,round((h.vbl - h.vblt0),3),h.trialNum,'.','Color',[0.6 0.6 0.6]);
                   end
                   start(h.tRefractory);
        end
    end
end

function refStart(~,~)
    global h
    % 安全停止lick counter，避免在已停止状态下再次停止
    if strcmp(h.tLickCounter.Running,'on')
        stop(h.tLickCounter);
    end
    %disp('timer 1 stopped by timer-2 StartFcn');
    writeDigitalPin(h.a,'D9',0);
end

function refEnd(~,~)
    global h
    % 安全启动lick counter，避免重复启动
    if ~strcmp(h.tLickCounter.Running,'on')
        start(h.tLickCounter);
    end
    % 安全停止refractory timer
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end
end

% GPU预热和缓存优化函数
function warmUpGPU()
    global h
    fprintf('正在预热GPU和优化缓存...\n');
    
    % 执行一系列渲染操作来预热GPU
    for i = 1:10
        % 绘制Gabor纹理进行预热
        Screen('DrawTexture', h.window, h.gratingtex, [], [], 0, [], [], [], [],...
              kPsychDontDoRotation, h.propertiesMat');
        Screen('Flip', h.window);
        
        % 更新相位以测试纹理动态更新
        h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
    end
    
    % 清屏并最后翻转
    Screen('FillRect', h.window, h.grey);
    Screen('Flip', h.window);
    
    % 重置相位
    h.propertiesMat(1) = 0;
    
    fprintf('GPU预热完成！\n');
end
