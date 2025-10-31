%% Enhanced Stage 3 for Two-Photon Imaging Synchronization
% Modified from stage3_Demo240119.m and stage2_RandITI240116.m
% 
% Key improvements for two-photon imaging:
% 1. 采用stage2的优化视觉刺激和UI设计
% 2. 保持stage3的target/non-target discrimination机制
% 3. 高精度时间戳记录，便于与双光子成像同步
% 4. 增强的数据输出格式，包含详细的事件时间轴
% 5. TTL信号输出功能，用于硬件同步
% 6. 实时lick rate热图监测
%
% Date: 2025-10-28
% Author: Enhanced for two-photon compatibility

timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor with enhanced performance settings
global h 
PsychDefaultSetup(2);

% 为双光子成像优化的Screen偏好设置 - 最高精度时间同步
Screen('Preference', 'ConserveVRAM', 4096); % 优化显存使用
Screen('Preference', 'VBLTimestampingMode', 4); % 最高精度VBL时间戳
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

% 预热显卡和优化缓存，同时清理任何残留的PTB缓冲区
try
    % 清理任何之前的PTB状态
    Screen('Close');  % 关闭所有纹理
    Screen('CloseAll'); % 确保清理
catch
    % 忽略清理错误
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
    warning('Steinberg UR12 device not found, using default audio device.');
    device_id = 1; % 使用默认音频设备
end

h.audioHandle = PsychPortAudio('Open',device_id-1, 1, 1, h.sampleF, 2);
PsychPortAudio('Volume', h.audioHandle, 0.1);      %auditory cue volume 

%pre-allocate audio buffer
[myBeep, samplingRate] = MakeBeep(10000, 0.1, h.sampleF);
buffer = [myBeep;myBeep];
PsychPortAudio('FillBuffer',h.audioHandle,buffer);

%% start communication with buffer management
try
    h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
    
    % 等待Arduino初始化完成
    pause(2);
    
    fprintf('Arduino connected successfully.\n');
catch ME
    error('Failed to connect to Arduino: %s', ME.message);
end

h.sensorPin = 'D13';
h.waterPumpPin = 'D9';
h.airPumpPin = 'D3';
% 新增TTL输出引脚用于双光子同步
h.ttlTrialStartPin = 'D2';  % Trial开始TTL信号
h.ttlVisualStimPin = 'D4';  % 视觉刺激TTL信号
h.ttlRewardPin = 'D5';      % 奖励TTL信号

% 初始化所有输出引脚为LOW状态
writeDigitalPin(h.a, h.waterPumpPin, 0);
writeDigitalPin(h.a, h.airPumpPin, 0);
writeDigitalPin(h.a, h.ttlTrialStartPin, 0);
writeDigitalPin(h.a, h.ttlVisualStimPin, 0);
writeDigitalPin(h.a, h.ttlRewardPin, 0);

% PTB显卡缓冲区管理计数器
h.frameCounter = 0;
h.bufferClearInterval = 50; % 每50帧清理一次PTB缓冲区

%% Gabor presetting - 采用stage2的优化设置
% 使用与visual_stimulation_gui相同的设置方式
[h.width, h.height] = Screen('WindowSize', h.window);
h.gaborDimPix = max(h.width, h.height); % 使用屏幕尺寸

% Target/Non-target orientation parameters (from stage3)
h.contrast = 1.0;
h.phase = 0;
h.orientationTarget = 0; %vertical - target
h.orientationNontarget = 90; %horizontal - non-target

% 使用与stage2相同的空间频率设置
h.numCycles = 7;
h.freq = h.numCycles / h.gaborDimPix;

% 使用stage2的优化参数
h.sigma = h.gaborDimPix;
h.aspectRatio = 1;
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;

% 创建Gabor纹理
h.gratingtex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

% 属性矩阵
h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];

% 优化的相位增量和帧设置
h.phasePerFrame = 5 * pi;
h.waitframes = 1;

% 预计算一些常用值以提高性能
h.halfIfi = 0.5 * h.ifi;
h.frameAdvance = (h.waitframes - 0.5) * h.ifi;

%% Enhanced program variables for two-photon synchronization
nTrial = 300;
nTarget = 150;   % 150 target trials out of 300 total
h.randSeq = randomSequence(nTrial, nTarget);  % randomize target/non-target sequence
h.oriSequence = h.randSeq * 90;    % create orientation sequence (0° target, 90° non-target)

h.licktrial = 1;
h.lickdata = {};
h.data1 = zeros(nTrial, 10); % 扩展数据矩阵以包含更多信息

% Trial counters for UI (from stage3)
h.totalLickTimes = 0;
h.hitTrialNumber = 0;
h.FATrialNumber = 0;
h.CRTrialNumber = 0;
h.missTrialNumber = 0;
h.resultFlag = [];
h.correctRate = 0;

% Enhanced timing variables for two-photon synchronization
h.sessionStartTime = GetSecs; % 会话开始的绝对时间
h.trialTimestamps = []; % 每个trial的详细时间戳
h.eventLog = {}; % 详细事件日志

% 高精度时间记录系统 - 替换mLatency
h.preciseTimings = struct(); % 精确的时间记录结构

% Lick rate monitoring (from stage2) - 热图版本
h.lickRateTimeWindow = 0.1; % 100ms时间窗
h.lickRateTimeStep = 0.05; % 50ms步长
h.trialLickTimes = []; % 当前trial中所有lick的时间戳
h.lickRateTimeCourse = {}; % 保存每个trial的lick count时间过程
h.timeAxis = 0:h.lickRateTimeStep:15; % 扩展时间轴至15秒以适应更长的试验

% 热图数据矩阵
h.lickRateMatrix = zeros(nTrial, length(h.timeAxis));
h.heatmapHandle = [];

% 创建随机ITI时间长度 4-6秒
h.ITIperiod = 4 + rand([1 nTrial])*1;

%% display UI and get mouse information
infoUI();

%% Enhanced UI combining stage2 and stage3 features
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;
f = figure('Name','Enhanced Two-Photon Behavioral Monitor','Position',screenSize,...
    'Color',[0.95 0.95 0.95]);

% Trial counters (optimized layout)
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

% Performance counters (from stage3)
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.25 0.95 0.06 0.03],'String','Hit','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',12);
h.hitTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.25 0.92 0.06 0.025],'FontSize',12);
h.hitCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.25 0.89 0.06 0.02],'String','','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',8);

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.33 0.95 0.06 0.03],'String','Miss','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',12);
h.missTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.33 0.92 0.06 0.025],'FontSize',12);
h.missCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.33 0.89 0.06 0.02],'String','','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',8);

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.41 0.95 0.06 0.03],'String','FA','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',12);
h.FATrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.41 0.92 0.06 0.025],'FontSize',12);
h.FACounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.41 0.89 0.06 0.02],'String','','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',8);

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.49 0.95 0.06 0.03],'String','CR','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',12);
h.CRTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.49 0.92 0.06 0.025],'FontSize',12);
h.CRCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.49 0.89 0.06 0.02],'String','','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',8);

% Target/Non-target indicators (from stage3)
h.targetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.93 0.1 0.05],'String','TARGET','BackgroundColor',[1 1 0],...
    'FontSize',14,'Visible','off');
h.nontargetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.82 0.1 0.05],'String','NON-TARGET','BackgroundColor',[0.5 0.5 1],...
    'FontSize',14,'Visible','off');

% Timing parameters display
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.57 0.95 0.06 0.03],'String','VS','BackgroundColor',[1 1 1],...
    'FontSize',12);
h.tempF =uicontrol(f,'Style','edit','String','1','Units','normalized',...
    'Position',[0.57 0.92 0.06 0.025],'FontSize',12); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.65 0.95 0.06 0.03],'String','RW','BackgroundColor',[1 1 1],...
    'FontSize',12);
h.spatF =uicontrol(f,'Style','edit','String','4','Units','normalized',...
    'Position',[0.65 0.92 0.06 0.025],'FontSize',12); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.73 0.95 0.06 0.03],'String','ITI','BackgroundColor',[1 1 1],...
    'FontSize',12);
h.duration =uicontrol(f,'Style','edit','String','4-6','Units','normalized',...
    'Position',[0.73 0.92 0.06 0.025],'FontSize',12); 

% 将三个坐标轴图放在中间一行，为底部摄像头预留空间
% Trial raster plot - 左侧
h.trialRaster = axes(f,'Position',[0.05 0.45 0.28 0.35],'FontSize',12);
title('Trial Raster');
xlabel('Time in trial (s)');
ylabel('Trial number');
xlim([0 15]);
ylim([0 300]);
set(h.trialRaster,'Ydir','reverse');
hold on

% Performance plot - 中间
h.ratePlot = axes(f,'Position',[0.37 0.45 0.28 0.35],'FontSize',12);
title('Performance Monitor');
xlabel('Trial number');     
xlim([1 300]);
ylim([0 1]);
ylabel('Correct rate');
hold on

% Enhanced lick rate heatmap - 右侧
h.lickRateMonitor = axes(f,'Position',[0.69 0.45 0.28 0.35],'FontSize',12);
title('Lick Count Heatmap (times/100ms)');
xlabel('Time in trial (s)');     
xlim([0 15]);
ylabel('Trial number');
ylim([0 300]);
set(h.lickRateMonitor,'YDir','reverse');
hold on

% Camera feeds - 重新布局到底部一行
h.backCam = webcam(1);
h.frontCam = webcam(2);

backRes = str2double(strsplit(h.backCam.Resolution,'x'));
frontRes = str2double(strsplit(h.frontCam.Resolution,'x'));

backAspect = backRes(1) / backRes(2);
frontAspect = frontRes(1) / frontRes(2);

% 调整摄像头窗口大小和位置，放在底部
backHeight = 0.25;
backWidth = backAspect * backHeight * 0.8; % 稍微缩小宽度
frontHeight = 0.25;
frontWidth = frontAspect * frontHeight * 0.8;

% 左侧摄像头 - Back Camera
h.backCamUI = axes(f, 'Position', [0.15 0.05 backWidth backHeight]);
h.im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', h.backCamUI);
preview(h.backCam, h.im);
text(h.backCamUI, 0.5, -0.05, 'Back Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% 右侧摄像头 - Front Camera  
h.frontCamUI = axes(f, 'Position', [0.55 0.05 frontWidth frontHeight]);
h.im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', h.frontCamUI);
preview(h.frontCam, h.im2);
text(h.frontCamUI, 0.5, -0.05, 'Front Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

%% define the enhanced timers
h.tLickCounter = timer('ExecutionMode', 'fixedRate', 'Period', 0.01,...
                      'TimerFcn',@(src,event)pinStatusChanged);
h.tRefractory = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.1,...
    'StartFcn',@(src,event)refStart,...
    'TimerFcn',@(src,event)refEnd);
h.tAirpuff = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.2,...
       'TimerFcn',@(src,event)airpuffEnd);

%% 3 minutes countdown 
countDown;

%% Enhanced main experiment loop with precise timing
totalTic = tic;
h.inOrOutRW = 0;
h.inRWCounter = 0;
h.outRWCounter = 0;
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
h.trialNum = [];
h.rwTime = [];
h.trialLickTimes = [];

% 初始化高精度时间记录矩阵
h.preciseTimings.trialData = zeros(nTrial, 12); % 扩展为12列以包含更多时间信息
createPreciseTimingDocumentation(); % 显示时间记录系统说明

fprintf('=== Enhanced Two-Photon Behavioral Session Started ===\n');
fprintf('Session start time: %f\n', h.sessionStartTime);

start(h.tLickCounter);
for trialNum = 1:nTrial
    h.trialNum = trialNum;
    
    % 安全停止refractory timer
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end
    
    % 安全管理lick counter
    if strcmp(h.tLickCounter.Running,'on')
        stop(h.tLickCounter);
    end
    
    h.lickInRWOneTrial = 0;
    h.outRWCounterSingle = 0;
    h.trialLickTimes = [];
    
    % Trial type determination
    h.visiOri = h.oriSequence(h.trialNum);
    h.targetFlag = (h.visiOri == h.orientationTarget);
    
    % Display trial type indicator
    if h.targetFlag
        set(h.targetBox,'Visible','on');
    else
        set(h.nontargetBox,'Visible','on');
    end
    
    % Record trial start with high precision timestamp
    h.trialGlobalTic = tic;
    trialAbsoluteStartTime = GetSecs;
    h.trialTimestamps(trialNum).absoluteStart = trialAbsoluteStartTime;
    h.trialTimestamps(trialNum).relativeStart = trialAbsoluteStartTime - h.sessionStartTime;
    
    % 记录精确时间信息
    h.preciseTimings.trialData(trialNum,1) = trialAbsoluteStartTime; % Trial start (absolute)
    h.preciseTimings.trialData(trialNum,2) = trialAbsoluteStartTime; % 将在trial结束时更新为end time
    
    % TTL signal for trial start (two-photon sync)
    writeDigitalPin(h.a, h.ttlTrialStartPin, 1);
    WaitSecs(0.001); % 1ms TTL pulse
    writeDigitalPin(h.a, h.ttlTrialStartPin, 0);
    
    logEvent(sprintf('Trial %d Start (Target: %d)', trialNum, h.targetFlag), GetSecs);
    
    disp('================== NEW TRIAL ==================')
    fprintf('Trial #%d | Type: %s | Orientation: %d°\n', trialNum, ...
        iif(h.targetFlag, 'TARGET', 'NON-TARGET'), h.visiOri);
    set(h.totalTrialNumUI,'String',num2str(trialNum));
   
    % 记录trial间隔时间（如果不是第一个trial）
    if trialNum > 1 
        h.preciseTimings.trialData(trialNum,12) = trialAbsoluteStartTime - h.preciseTimings.trialData(trialNum-1,2);
    end

    % Phase 1: Trial Cue - 高精度时间记录
    cueStartTime = GetSecs;
    trialCue;
    cueEndTime = GetSecs;
    h.preciseTimings.trialData(trialNum,3) = cueStartTime; % Cue start (absolute)
    h.preciseTimings.trialData(trialNum,4) = cueEndTime - cueStartTime; % Cue duration

    % Phase 2: Post-cue Period
    postCueStartTime = GetSecs;
    postCuePeriod;
    postCueEndTime = GetSecs;
    h.preciseTimings.trialData(trialNum,5) = postCueEndTime - postCueStartTime; % Post-cue duration

    % Phase 3: Visual Stimulation
    vsStartTime = GetSecs;
    visiStim;
    vsEndTime = GetSecs;
    h.preciseTimings.trialData(trialNum,6) = vsStartTime; % VS start (absolute)
    h.preciseTimings.trialData(trialNum,7) = vsEndTime - vsStartTime; % VS duration

    % Phase 4: Response Window
    rwStartTime = GetSecs;
    responseWindow;
    rwEndTime = GetSecs;
    h.preciseTimings.trialData(trialNum,8) = rwStartTime; % RW start (absolute)
    h.preciseTimings.trialData(trialNum,9) = rwEndTime - rwStartTime; % RW duration

    % Phase 5: ITI Period
    itiStartTime = GetSecs;
    ITIperiod;
    itiEndTime = GetSecs;
    h.preciseTimings.trialData(trialNum,10) = itiEndTime - itiStartTime; % ITI duration
    
    % Calculate and plot lick rate for current trial
    lickAnalysisStartTime = GetSecs;
    calculateAndPlotLickRate(trialNum);
    h.preciseTimings.trialData(trialNum,11) = GetSecs - lickAnalysisStartTime; % Analysis time
    
    % Record trial end timestamp
    trialAbsoluteEndTime = GetSecs;
    h.trialTimestamps(trialNum).absoluteEnd = trialAbsoluteEndTime;
    h.trialTimestamps(trialNum).duration = toc(h.trialGlobalTic);
    
    % 更新精确时间记录
    h.preciseTimings.trialData(trialNum,2) = trialAbsoluteEndTime; % Trial end (absolute)
    
    % Performance calculation and data recording
    h.correctRate = (h.hitTrialNumber + h.CRTrialNumber) / trialNum;
    
    % Enhanced data recording
    h.data1(trialNum,1) = trialNum;
    h.data1(trialNum,2) = h.resultFlag;   % 1:Hit, 2:Miss, 3:FA, 4:CR
    h.data1(trialNum,3) = h.targetFlag;   % 1:Target, 0:Non-target
    h.data1(trialNum,4) = h.lickInRWOneTrial;
    h.data1(trialNum,5) = h.outRWCounterSingle;
    h.data1(trialNum,6) = h.trialTimestamps(trialNum).duration;
    h.data1(trialNum,7) = h.correctRate;
    h.data1(trialNum,8) = length(h.trialLickTimes); % Total licks in trial
    h.data1(trialNum,9) = h.trialTimestamps(trialNum).absoluteStart;
    h.data1(trialNum,10) = h.trialTimestamps(trialNum).absoluteEnd;
    
    % Update performance plot
    plot(h.ratePlot, trialNum, h.correctRate, '-ok');
    
    fprintf('Trial Summary: RW Licks=%d | Total Licks=%d | Performance=%.2f\n', ...
        h.lickInRWOneTrial, length(h.trialLickTimes), h.correctRate);
    
    % Reset trial type indicators
    set(h.targetBox,'Visible','off');
    set(h.nontargetBox,'Visible','off');
    
    % Reset UI background colors
    resetUIColors();
    
    % PTB缓冲区管理（每10个试次）
    if mod(trialNum, 10) == 0
        try
            % 清理PTB纹理和缓冲区
            Screen('Close'); % 关闭未使用的纹理
            fprintf('PTB buffers cleaned at trial %d\n', trialNum);
        catch
            fprintf('Warning: Failed to clean PTB buffers at trial %d\n', trialNum);
        end
    end
end

%% Session cleanup and data saving
stop(h.tLickCounter);
delete(h.tLickCounter);
delete(h.tRefractory);
delete(h.tAirpuff);

disp('=================== SESSION COMPLETED ===================')
totalTime = toc(totalTic);

% Calculate total statistics
totalLickTimes = sum([h.data1(:,8)]);
sessionEndTime = GetSecs;

PsychPortAudio('Close',h.audioHandle);

fprintf('\n=== SESSION STATISTICS ===\n');
fprintf('Total time: %.1f minutes (%.1f seconds)\n', totalTime/60, totalTime);
fprintf('Total trials: %d\n', trialNum);
fprintf('Total licks: %d\n', totalLickTimes);
fprintf('Hit trials: %d\n', h.hitTrialNumber);
fprintf('Miss trials: %d\n', h.missTrialNumber);
fprintf('FA trials: %d\n', h.FATrialNumber);
fprintf('CR trials: %d\n', h.CRTrialNumber);
fprintf('Final correct rate: %.2f%%\n', h.correctRate * 100);

% Enhanced data structure for two-photon analysis
twoPhotonData.sessionInfo.startTime = h.sessionStartTime;
twoPhotonData.sessionInfo.endTime = sessionEndTime;
twoPhotonData.sessionInfo.duration = totalTime;
twoPhotonData.sessionInfo.mouseID = h.mouseID;
twoPhotonData.sessionInfo.ifi = h.ifi; % Screen refresh rate for synchronization

twoPhotonData.trialData = h.data1;
twoPhotonData.trialTimestamps = h.trialTimestamps;
twoPhotonData.eventLog = h.eventLog;
twoPhotonData.lickData = h.lickdata;
twoPhotonData.lickRateTimeCourse = h.lickRateTimeCourse;
twoPhotonData.lickRateMatrix = h.lickRateMatrix;
twoPhotonData.orientationSequence = h.oriSequence;
twoPhotonData.ITIperiods = h.ITIperiod;
twoPhotonData.preciseTimings = h.preciseTimings; % 替换mLatency为精确时间记录

% Save with timestamp and mouse ID
saveFilename = sprintf('TwoPhoton_%s_%s_stage3.mat', ...
    datestr(now, 'yyyymmdd_HHMMSS'), char(h.mouseID{1}));
save(saveFilename, 'twoPhotonData', 'h');

fprintf('Data saved to: %s\n', saveFilename);

% Export synchronization data for two-photon analysis
exportSyncData(saveFilename);

sca;

%% Enhanced Functions

function calculateAndPlotLickRate(trialNum)
global h
    % Calculate lick count in sliding windows
    lickCount = zeros(size(h.timeAxis));
    
    for i = 1:length(h.timeAxis)
        windowStart = h.timeAxis(i) - h.lickRateTimeWindow/2;
        windowEnd = h.timeAxis(i) + h.lickRateTimeWindow/2;
        
        licksInWindow = sum(h.trialLickTimes >= windowStart & h.trialLickTimes <= windowEnd);
        lickCount(i) = licksInWindow;
    end
    
    h.lickRateTimeCourse{trialNum} = lickCount;
    h.lickRateMatrix(trialNum, :) = lickCount;
    
    updateLickRateHeatmap(trialNum);
end

function updateLickRateHeatmap(trialNum)
global h
    dataToShow = h.lickRateMatrix(1:trialNum, :);
    
    if trialNum == 1
        h.heatmapHandle = imagesc(h.lickRateMonitor, h.timeAxis, 1:trialNum, dataToShow);
        colormap(h.lickRateMonitor, 'hot');
        colorbar(h.lickRateMonitor);
        
        set(h.lickRateMonitor, 'YDir', 'reverse');
        xlabel(h.lickRateMonitor, 'Time in trial (s)');
        ylabel(h.lickRateMonitor, 'Trial number');
        title(h.lickRateMonitor, 'Lick Count Heatmap (times/100ms)');
        
        xlim(h.lickRateMonitor, [0 15]);
        ylim(h.lickRateMonitor, [0.5 300.5]);
        
        % Add phase markers
        hold(h.lickRateMonitor, 'on');
        plot(h.lickRateMonitor, [1 1], [0.5 300.5], 'w--', 'LineWidth', 1); % Post-cue end
        plot(h.lickRateMonitor, [2 2], [0.5 300.5], 'w--', 'LineWidth', 1); % VS end
        plot(h.lickRateMonitor, [6 6], [0.5 300.5], 'w--', 'LineWidth', 1); % RW end
        text(h.lickRateMonitor, 0.5, 10, 'Cue', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 1.5, 10, 'VS', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 4, 10, 'RW', 'Color', 'white', 'FontSize', 8);
        text(h.lickRateMonitor, 10, 10, 'ITI', 'Color', 'white', 'FontSize', 8);
        hold(h.lickRateMonitor, 'off');
    else
        set(h.heatmapHandle, 'CData', dataToShow, 'YData', 1:trialNum);
    end
    
    caxis(h.lickRateMonitor, [0 3]);
    drawnow;
end

function logEvent(eventString, timestamp)
global h
    eventIndex = length(h.eventLog) + 1;
    h.eventLog{eventIndex}.event = eventString;
    h.eventLog{eventIndex}.absoluteTime = timestamp;
    h.eventLog{eventIndex}.relativeTime = timestamp - h.sessionStartTime;
    
    fprintf('EVENT: %s @ %.3f s\n', eventString, timestamp - h.sessionStartTime);
end

function resetUIColors()
global h
    set(h.hitCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
    set(h.missCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
    set(h.FACounterUI,'BackgroundColor',[0.9 0.9 0.9]);
    set(h.CRCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
end

function exportSyncData(filename)
global h
    % Export synchronization file for two-photon analysis software
    [filepath, name, ~] = fileparts(filename);
    syncFilename = fullfile(filepath, [name '_sync.csv']);
    
    % Create sync table with key events
    syncData = [];
    for i = 1:length(h.trialTimestamps)
        if ~isempty(h.trialTimestamps(i).absoluteStart)
            syncData = [syncData; i, h.trialTimestamps(i).absoluteStart, ...
                       h.trialTimestamps(i).absoluteEnd, h.trialTimestamps(i).duration];
        end
    end
    
    % Write to CSV
    if ~isempty(syncData)
        T = table(syncData(:,1), syncData(:,2), syncData(:,3), syncData(:,4), ...
            'VariableNames', {'TrialNumber', 'StartTime', 'EndTime', 'Duration'});
        writetable(T, syncFilename);
        fprintf('Synchronization data exported to: %s\n', syncFilename);
    end
end

function result = iif(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end

function infoUI(~,~)
global h
    prompt = {'Mouse ID', 'Training Stage', 'Day Number', 'Save Directory', 'Experimenter'};
    dlgtitle = 'Two-Photon Session Information';
    dims = [1 50];
    definput = {'', 'Stage3_TwoPhoton', '', pwd, ''};
    h.mouseID = inputdlg(prompt, dlgtitle, dims, definput);
end

function countDown(~,~)
    countdownDuration = 180; % 3 minutes
    disp('========== Two-Photon Session Countdown ==========')
    disp('Check: Mouse position, lick spout, two-photon alignment!')
    
    for remainingSeconds = countdownDuration:-1:0
        fprintf('Session starts in: %d seconds\n', remainingSeconds);
        pause(1);
    end
    disp('========== Two-Photon Session Started! ==========')
end

%% Trial procedure functions

function trialCue(~,~)
global h
    % Record cue start time
    cueStartTime = GetSecs;
    logEvent('Auditory Cue Start', cueStartTime);
    
    PsychPortAudio('Start',h.audioHandle,1,0,1);
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop',h.audioHandle,1,1);
    
    cueEndTime = GetSecs;
    logEvent('Auditory Cue End', cueEndTime);
    
    disp('>> Trial Cue completed! Post-cue period starts!')
end

function postCuePeriod(~,~)
global h
    h.inOrOutRW = -1;  % Post-cue period marker
    
    postCueStartTime = GetSecs;
    logEvent('Post-Cue Period Start', postCueStartTime);
    
    % Safely start lick counter
    if ~strcmp(h.tLickCounter.Running,'on')
        start(h.tLickCounter); 
    end
    
    h.postCueTime = tic;
    while toc(h.postCueTime) <= 1  % Fixed 1 second
        WaitSecs(0.001); % Prevent busy waiting
    end
    
    postCueEndTime = GetSecs;
    logEvent('Post-Cue Period End', postCueEndTime);
    
    disp('>> Post-cue period completed! Visual stimulation starts!')
end

function visiStim(~,~)
global h
    h.inOrOutRW = 0; % Visual stimulation period marker
    
    try
        % TTL signal for visual stimulus start (two-photon sync)
        writeDigitalPin(h.a, h.ttlVisualStimPin, 1);
        
        visStimStartTime = GetSecs;
        logEvent(sprintf('Visual Stim Start (Ori: %d°)', h.visiOri), visStimStartTime);
        
        % PTB缓冲区管理 - 定期清理显卡缓冲区
        h.frameCounter = h.frameCounter + 1;
        if h.frameCounter > h.bufferClearInterval
            % 清理PTB内部缓冲区
            Screen('Close'); % 关闭不必要的纹理
            % 强制垃圾回收
            Screen('Preference', 'ConserveVRAM', 4096);
            h.frameCounter = 0;
        end
        
        % Enhanced visual stimulation with precise timing and buffer management
        Screen('FillRect', h.window, h.grey);
        
        % 同步显卡状态
        Screen('DrawingFinished', h.window); % 确保之前的绘图完成
        
        vbl = Screen('Flip', h.window);
        vblt0 = vbl;
        startTime = vbl;
        
        h.propertiesMat(1) = 0; % Reset phase
        
        VSLength = 1; % 1 second visual stimulation
        frameCount = 0;
        maxFrames = ceil(VSLength / h.ifi); % 预计算最大帧数
        
        while (vbl - vblt0) <= VSLength && frameCount < maxFrames
            frameCount = frameCount + 1;
            
            % 安全的纹理绘制
            try
                % Draw Gabor with current orientation
                Screen('DrawTexture', h.window, h.gratingtex, [], [], h.visiOri, [], [], [], [],...
                      kPsychDontDoRotation, h.propertiesMat');
                
                % 确保绘图完成
                Screen('DrawingFinished', h.window);
                
                % Precise timing for next flip
                nextFlipTime = startTime + frameCount * h.ifi;
                vbl = Screen('Flip', h.window, nextFlipTime - 0.5 * h.ifi);
                
                % Update phase for motion
                h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
                
            catch flipError
                % 如果出现缓冲区错误，尝试清理并重试
                if contains(flipError.message, 'buffer') || contains(flipError.message, 'full')
                    fprintf('PTB buffer issue detected, clearing buffers...\n');
                    Screen('DrawingFinished', h.window, 1); % 强制完成所有绘图
                    WaitSecs(0.001); % 短暂等待
                    continue;
                else
                    rethrow(flipError);
                end
            end
        end
        
        % Clear screen and final flip
        Screen('FillRect', h.window, h.grey);
        Screen('DrawingFinished', h.window);
        Screen('Flip', h.window);
        
        % TTL signal for visual stimulus end
        writeDigitalPin(h.a, h.ttlVisualStimPin, 0);
        
        visStimEndTime = GetSecs;
        logEvent('Visual Stim End', visStimEndTime);
        
        h.inOrOutRW = 1; % Now in response window
        disp('>> Visual stimulation completed! Response window starts!')
        
    catch ME
        % 错误处理
        fprintf('Error in visual stimulation: %s\n', ME.message);
        
        % 尝试恢复
        try
            Screen('FillRect', h.window, h.grey);
            Screen('Flip', h.window);
            writeDigitalPin(h.a, h.ttlVisualStimPin, 0);
        catch
            % 如果恢复失败，记录错误但继续
            fprintf('Failed to recover from visual stimulation error\n');
        end
        
        h.inOrOutRW = 1; % 确保状态正确
        rethrow(ME);
    end
end

function responseWindow(~,~)
global h
    rwStartTime = GetSecs;
    logEvent('Response Window Start', rwStartTime);
    
    h.rwTime = tic;
    h.rwLimit = 4; % 4 second response window
    
    while toc(h.rwTime) <= h.rwLimit
        WaitSecs(0.001); % Prevent busy waiting
    end
    
    % Update lick trial counter if there were licks in RW
    if h.lickInRWOneTrial > 0
        h.licktrial = h.licktrial + 1;
    end
    
    h.inOrOutRW = 0; % Reset marker
    
    rwEndTime = GetSecs;
    logEvent('Response Window End', rwEndTime);
    
    disp('>> Response window completed! ITI starts!')
end

function ITIperiod(~,~)
global h
    h.inOrOutRW = 2; % ITI period marker
    
    itiStartTime = GetSecs;
    logEvent('ITI Start', itiStartTime);
    
    % Determine trial outcome if no licks in RW
    if h.lickInRWOneTrial == 0
        if h.targetFlag    % Miss
            h.missTrialNumber = h.missTrialNumber + 1;
            h.resultFlag = 2;
            set(h.missTrialNumUI,'String',num2str(h.missTrialNumber));
            set(h.missCounterUI,'BackgroundColor',[1 0 0]);
            logEvent('Trial Outcome: Miss', GetSecs);
        else   % Correct Rejection
            h.CRTrialNumber = h.CRTrialNumber + 1;
            h.resultFlag = 4;
            set(h.CRTrialNumUI,'String',num2str(h.CRTrialNumber));
            set(h.CRCounterUI,'BackgroundColor',[0 1 0]);
            logEvent('Trial Outcome: Correct Rejection', GetSecs);
        end
    end
    
    h.ITITime = tic;
    itiDuration = h.ITIperiod(h.trialNum);
    
    while toc(h.ITITime) <= itiDuration
        WaitSecs(0.01); % Reduced CPU usage
    end
    
    itiEndTime = GetSecs;
    logEvent(sprintf('ITI End (Duration: %.1fs)', itiDuration), itiEndTime);
    
    fprintf('>> ITI completed! Duration: %.1f seconds\n', itiDuration);
end

function pinStatusChanged(~,~)
global h
    % Enhanced lick detection with precise timing and buffer management
    if ~strcmp(h.tLickCounter.Running,'on')
        return; % Safety check
    end
    
    % 移除错误的Arduino缓冲区管理代码
    % PTB缓冲区问题已在visiStim函数中处理
    
    RWflag = h.inOrOutRW;
    lickFlag = true;
    trialFlag = true;
    lickTimesReporter = 0;
    
    try
        pinValue = readDigitalPin(h.a, h.sensorPin);
    catch ME
        % 如果读取失败，记录错误但继续
        fprintf('Warning: Arduino read failed: %s. Skipping this read.\n', ME.message);
        return;
    end
    
    if pinValue == true && lickFlag == true
        lickTime = GetSecs;
        trialElapsedTime = toc(h.trialGlobalTic);
        h.trialLickTimes = [h.trialLickTimes, trialElapsedTime];
        
        switch RWflag
            case 1 % Response Window
                h.inRWCounter = h.inRWCounter + 1;
                
                if h.targetFlag  % Target trial
                    if h.lickInRWOneTrial < 1  % First lick = HIT
                        lickTimeInRW = toc(h.rwTime);
                        if lickTimeInRW > 10
                            plot(h.trialRaster, 1.0005, h.trialNum, '.g');  
                        else
                            plot(h.trialRaster, lickTimeInRW + 1, h.trialNum, '.g'); 
                        end
                        
                        % Reward delivery with TTL
                        writeDigitalPin(h.a, h.waterPumpPin, 1);
                        writeDigitalPin(h.a, h.ttlRewardPin, 1);
                        WaitSecs(0.001);
                        writeDigitalPin(h.a, h.ttlRewardPin, 0);
                        
                        h.hitTrialNumber = h.hitTrialNumber + 1; 
                        h.resultFlag = 1;
                        set(h.hitTrialNumUI,'String',num2str(h.hitTrialNumber));
                        set(h.hitCounterUI,'BackgroundColor',[0 1 0]);
                        
                        logEvent(sprintf('HIT (RT: %.3fs)', lickTimeInRW), lickTime);
                    else
                        % Additional licks in RW
                        plot(h.trialRaster, toc(h.rwTime) + 1, h.trialNum, '.', 'Color', [0.6 0.6 0.6]); 
                    end
                    
                else  % Non-target trial
                    if h.lickInRWOneTrial < 1  % First lick = FALSE ALARM
                        lickTimeInRW = toc(h.rwTime);
                        if lickTimeInRW > 10
                            plot(h.trialRaster, 1.001, h.trialNum, '.r');  
                        else
                            plot(h.trialRaster, lickTimeInRW + 1, h.trialNum, '.r');
                        end
                        
                        % Air puff punishment
                        writeDigitalPin(h.a, h.airPumpPin, 1);
                        h.FATrialNumber = h.FATrialNumber + 1;
                        h.resultFlag = 3;
                        set(h.FATrialNumUI,'String',num2str(h.FATrialNumber));   
                        set(h.FACounterUI,'BackgroundColor',[1 0 0]);
                        start(h.tAirpuff);
                        
                        logEvent(sprintf('FALSE ALARM (RT: %.3fs)', lickTimeInRW), lickTime);
                        disp('FALSE ALARM! Air puff delivered.');
                    else
                        % Additional licks in RW
                        plot(h.trialRaster, toc(h.rwTime) + 1, h.trialNum, '.', 'Color', [0.6 0.6 0.6]); 
                    end
                end
                
                start(h.tRefractory);
                lickFlag = false;
                
                % Record lick data
                if lickFlag == false && trialFlag == true
                    trialFlag = false;
                    fprintf('Lick trial = %d ', h.licktrial);
                    set(h.lickTrialNumUI,'String',num2str(h.licktrial));
                    
                    lickTimesReporter = lickTimesReporter + 1;
                    h.lickInRWOneTrial = h.lickInRWOneTrial + lickTimesReporter;
                    tocReporter = toc(h.rwTime);
                    
                    fprintf('Lick #%d @ %.3fs in RW\n', h.lickInRWOneTrial, tocReporter);
                    h.lickdata{h.licktrial}(h.lickInRWOneTrial,1) = h.lickInRWOneTrial;
                    h.lickdata{h.licktrial}(h.lickInRWOneTrial,2) = tocReporter;
                end
                
            case -1 % Post-cue period - early lick resets timer
                h.postCueTime = tic;
                logEvent('Early Lick - Reset Post-Cue', lickTime);
                disp('!! Early lick! Reset post-cue period timer!');
                start(h.tRefractory);
                
            case 0 % Visual stimulation period
                if exist('h.vbl','var') && exist('h.vblt0','var')
                    if (h.vbl - h.vblt0) <= 1
                        plot(h.trialRaster, (h.vbl - h.vblt0), h.trialNum, '.', 'Color', [0.6 0.6 0.6]);
                        logEvent('Lick during Visual Stim', lickTime);
                    end
                end
                start(h.tRefractory);
                
            case 2 % ITI period
                plot(h.trialRaster, toc(h.ITITime) + 1 + h.rwLimit, h.trialNum, '.', 'Color', [0.6 0.6 0.6]);
                h.outRWCounter = h.outRWCounter + 1;
                h.outRWCounterSingle = h.outRWCounterSingle + 1;
                logEvent('ITI Lick', lickTime);
                start(h.tRefractory);
        end
    end
end

function refStart(~,~)
global h
    if strcmp(h.tLickCounter.Running,'on')
        stop(h.tLickCounter);
    end
    writeDigitalPin(h.a, h.waterPumpPin, 0);
end

function refEnd(~,~)
global h
    if ~strcmp(h.tLickCounter.Running,'on')
        start(h.tLickCounter);
    end
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end
end

function airpuffEnd(~,~)
global h
    writeDigitalPin(h.a, h.airPumpPin, 0);
    h.rwLimit = toc(h.rwTime) + 7; % Extend timeout
    stop(h.tAirpuff);
end

function createPreciseTimingDocumentation()
    % 精确时间记录系统说明
    % h.preciseTimings.trialData 矩阵列说明：
    % 第1列: Trial开始时间 (绝对时间戳)
    % 第2列: Trial结束时间 (绝对时间戳)  
    % 第3列: Cue开始时间 (绝对时间戳)
    % 第4列: Cue持续时间 (秒)
    % 第5列: Post-cue持续时间 (秒)
    % 第6列: 视觉刺激开始时间 (绝对时间戳)
    % 第7列: 视觉刺激持续时间 (秒)
    % 第8列: 反应窗口开始时间 (绝对时间戳)
    % 第9列: 反应窗口持续时间 (秒)
    % 第10列: ITI持续时间 (秒)
    % 第11列: Lick分析处理时间 (秒)
    % 第12列: Trial间隔时间 (秒)
    
    fprintf('Precise timing system initialized with 12-column data matrix.\n');
    fprintf('All absolute timestamps are in GetSecs format for microsecond precision.\n');
end

function seq = randomSequence(n, m)
    % Generate randomized sequence with no more than 3 consecutive repeats
    seq = zeros(n, 1);
    seq(1:m) = 1;
    seq(randperm(n)) = seq;
    
    nPermit = 3;
    iSame = 0;
    for i = 1:n-1
        if seq(i) == seq(i + 1)
            iSame = iSame + 1;
        else
            iSame = 0;
        end
        
        j = 0;
        while iSame > nPermit - 1
            temp = seq(end - j); 
            seq(end - j) = seq(i + 1); 
            seq(i + 1) = temp;
            
            if seq(i) == seq(i + 1)
                j = j + 1;
            else
                j = 0;
                iSame = 0;
            end
        end
    end
end