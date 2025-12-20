%% modified program from stageTwo.m - 整合stage2 UI设计和stage3 conditioning机制
%
timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor, h the gray color background
global h 
PsychDefaultSetup(2);

% 添加性能优化的Screen偏好设置（来自stage2）
Screen('Preference', 'ConserveVRAM', 4096); 
Screen('Preference', 'VBLTimestampingMode', 4); 
Screen('Preference', 'SkipSyncTests', 0); 
Screen('Preference', 'VisualDebugLevel', 0); 
Screen('Preference', 'SuppressAllWarnings', 1); 

Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);
h.screenNumber = max(Screen('Screens'));
h.white = WhiteIndex(h.screenNumber);
h.grey = h.white / 2;

% 优化窗口打开参数（来自stage2）
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution); 
h.ifi = Screen('GetFlipInterval',h.window); 

% 设置最高优先级和实时调度
h.topPriorityLevel = MaxPriority(h.window); 
Priority(h.topPriorityLevel);

% 预热显卡和优化缓存（来自stage2）
for i = 1:5
    Screen('Flip', h.window);
end

%% initialize sound configuration
InitializePsychSound;

%open psych-audio port - 使用stage2的设备查找方式
h.sampleF = 48000;
deviceList = PsychPortAudio('GetDevices');

% find the device ID for the Steinberg UR12（来自stage2）
device_id = [];
for i = 1:length(deviceList)
    if contains(deviceList(i).DeviceName, 'UR12') || contains(deviceList(i).DeviceName, 'Steinberg')
        device_id = i;
        disp(['Audio device found with Device ID: ' num2str(device_id)]);
        break;
    end
end
if isempty(device_id)
    % 回退方案：使用原stage3的设备ID
    %device_id = 13; % 原stage3使用的设备ID
    warning('Steinberg UR12 device not found. Using fallback device ID: %d', device_id);
end

h.audioHandle = PsychPortAudio('Open', device_id-1, 1, 1, h.sampleF, 2);
PsychPortAudio('Volume', h.audioHandle, 0.02); % 保持stage3的音量设置

%pre-allocate audio buffer
[myBeep, samplingRate] = MakeBeep(10000, 0.1, h.sampleF);
buffer = [myBeep;myBeep];
PsychPortAudio('FillBuffer',h.audioHandle,buffer);

%% start communication
h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
h.sensorPin = 'D13';
h.waterPumpPin = 'D9';
h.airPumpPin = 'D3'; % 保留stage3的air pump

%% Gabor presetting - 使用stage2的优化设置但保留stage3的参数
[h.width, h.height] = Screen('WindowSize', h.window);
h.gaborDimPix = h.windowRect(4)*2; % 保留stage3的尺寸设置

%center of display position（保留stage3设置）
h.center = [(h.width-h.height)/2,0,(h.width+h.height)/2,h.height];

%other parameters - 保留stage3的参数
h.sigma = h.gaborDimPix;
h.oriTarget = 0; % vertical（保留stage3设置）
h.oriNonTarget = 90; % horizontal（保留stage3设置）
h.contrast = 1;
h.aspectRatio = 1;
h.phase = 0; 

%spatial frequency - 保留stage3的设置
h.numCycles = 7;
h.freq = h.numCycles / h.gaborDimPix;

%make procedural gabor texture - 使用stage2的优化方式
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;
h.gabortex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

%make the property matrix
h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];
updateVbl;
h.waitframes = 1;
h.phasePerFrame = 5 * pi;  

% 预计算常用值以提高性能（来自stage2）
h.halfIfi = 0.5 * h.ifi;
h.frameAdvance = (h.waitframes - 0.5) * h.ifi;

%% program variables - 保留stage3的conditioning机制
nTrial = 300;
nTarget = 150;   
h.randSeq = randomSequence(nTrial,nTarget);  
h.oriSequence = h.randSeq * 90;    
h.licktrial = 1;
h.lickdata = {};
h.data1 = zeros(nTrial,8); % 扩展为8列以存储更多信息

%trial counters for UI - 保留stage3的所有计数器
h.totalLickTimes = 0;
h.hitTrialNumber = 0;
h.FATrialNumber = 0;
h.CRTrialNumber = 0;
h.missTrialNumber = 0;
h.resultFlag = [];
h.correctRate = 0;

% Lick rate monitoring variables - 引入stage2的热图功能
h.lickRateTimeWindow = 0.1; 
h.lickRateTimeStep = 0.05; 
h.trialLickTimes = []; 
h.lickRateTimeCourse = {}; 
h.timeAxis = 0:h.lickRateTimeStep:15; % 扩展到15秒以适应stage3的更长trial

% 热图数据矩阵
h.lickRateMatrix = zeros(nTrial, length(h.timeAxis));
h.heatmapHandle = [];

%create random ITI time length ranging from 4 ~ 6s.
h.ITIperiod = 4 + rand([1 nTrial])*2;

%% UI, waiting for the initialization
infoUI();

%% UI - 使用stage2的优化UI布局但保留stage3的所有功能
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);          

% 调整状态栏布局 - 更紧凑的单行设计
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.02 0.96 0.06 0.025],'String','Total Trial','BackgroundColor',[0 1 1],...
    'FontSize',9);
 h.totalTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.08 0.96 0.04 0.02],'FontSize',9); 
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.13 0.96 0.06 0.025],'String','Lick Trial','BackgroundColor',[0 1 1],...
    'FontSize',9);
 h.lickTrialNumUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.19 0.96 0.04 0.02],'FontSize',9);

 % 结果计数器 - 调整为单行紧凑布局
h.hitCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.25 0.96 0.08 0.025],'String','Hit Trial','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',9);
 h.hitTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.33 0.96 0.04 0.02],'FontSize',9); 

 h.missCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.38 0.96 0.08 0.025],'String','Miss Trial','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',9);
 h.missTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.46 0.96 0.04 0.02],'FontSize',9); 

 h.FACounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.51 0.96 0.08 0.025],'String','FA Trial','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',9);
 h.FATrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.59 0.96 0.04 0.02],'FontSize',9); 

 h.CRCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.64 0.96 0.08 0.025],'String','CR Trial','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',9);
 h.CRTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.72 0.96 0.04 0.02],'FontSize',9); 

 %target and nontarget box indicator - 调整到右上角单行
 h.targetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.78 0.96 0.06 0.025],'String','Target','BackgroundColor',[1 1 0],...
    'FontSize',9,'Visible','off');
 h.nontargetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.85 0.96 0.08 0.025],'String','Nontarget','BackgroundColor',[1 1 0],...
    'FontSize',9,'Visible','off');

 % 移到第二行，为状态栏腾出空间
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.02 0.93 0.03 0.025],'String','VS','BackgroundColor',[1 1 1],...
     'FontSize',9);
 h.tempF =uicontrol(f,'Style','edit','String','1','Units','normalized',...
    'Position',[0.05 0.93 0.03 0.02],'FontSize',9); 
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.09 0.93 0.03 0.025],'String','RW','BackgroundColor',[1 1 1],...
     'FontSize',9);
 h.spatF =uicontrol(f,'Style','edit','String','4','Units','normalized',...
    'Position',[0.12 0.93 0.03 0.02],'FontSize',9); 

 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.16 0.93 0.03 0.025],'String','ITI','BackgroundColor',[1 1 1],...
     'FontSize',9);
 h.duration =uicontrol(f,'Style','edit','String','4~6','Units','normalized',...
    'Position',[0.19 0.93 0.04 0.02],'FontSize',9); 

% trial raster - 调整位置以适应新的状态栏布局
h.trialRaster = axes(f,'Position',[0.05 0.45 0.28 0.4],'FontSize',12);
title('Trial Raster');
xlabel('seconds');
ylabel('trial number');
xlim([0 15]); % 适应stage3的更长时间
ylim([1 300]);
set(h.trialRaster,'Ydir','reverse');
hold on

% success rate monitor - 扩大到原来热图的位置
h.ratePlot = axes(f,'Position',[0.35 0.5 0.6 0.35],'FontSize',12);
title('Correct Rate Plot');
xlabel('trial number');     
xlim([1 300]);
ylim([0 1]);
ylabel('correct rate');
hold on

% lick rate monitor - 保留数据收集但不显示热图
% h.lickRateMonitor = axes(f,'Position',[0.65 0.5 0.3 0.35],'FontSize',12);
% title('Lick Count Heatmap (times/trial/100ms)');
% xlabel('Time in trial (s)');     
% xlim([0 15]); % 适应stage3的时间范围
% ylabel('Trial number');
% ylim([0 300]); 
% set(h.lickRateMonitor,'YDir','reverse'); 
% hold on

%open cameras - 使用stage2的双摄像头功能
try
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
catch
    warning('Camera initialization failed. Continuing without cameras.');
end

%% define the timers
h.tLickCounter = timer('ExecutionMode', 'fixedRate', 'Period', 0.01,...
                             'TimerFcn',@(src,event)pinStatusChanged);
%global lick report timer with 10ms period
h.tRefractory = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.1,...
    'StartFcn',@(src,event)refStart,...
    'TimerFcn',@(src,event)refEnd);
h.tAirpuff = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.2,...
       'TimerFcn',@(src,event)airpuffEnd);

% 新增：水泵控制专用timer，确保奖励持续时间
h.tWaterPump = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.15,...
       'TimerFcn',@(src,event)waterPumpEnd);

%% 3 minutes countdown 
countDown;

%% start tic - 优化时间戳记录系统
totalTic = tic;
h.inOrOutRW = [];
h.inRWCounter = 0;
h.outRWCounter = 0;
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
h.earlyTrialNumber = 0; 
h.trialNum = [];
h.rwTime = [];



mLatency = zeros(nTrial,10);   % 扩展为10列存储更多时间信息

%% main loop
start(h.tLickCounter);   
for trialNum = 1:nTrial
    h.trialNum = trialNum;
    
    % 安全停止timers
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end
    if strcmp(h.tLickCounter.Running,'on')
        stop(h.tLickCounter);
    end

    % 重置trial变量
    h.outRWCounterSingle = 0;
    h.lickInRWOneTrial = 0;
    h.trialLickTimes = []; 
    
    % 获取target/non-target信息
    h.visiOri = h.oriSequence(h.trialNum); 
    h.targetFlag = (h.visiOri == h.oriTarget); 
    
    % 显示target/non-target指示器
    if h.targetFlag     
       set(h.targetBox,'Visible','on');
    elseif ~h.targetFlag
        set(h.nontargetBox,'Visible','on');
    end

    % 记录trial开始的高精度时间戳
    h.trialGlobalTic = tic;
    
    disp('------------------new trial-----------------------')
    fprintf('trialNum = %s, Target = %d\n', num2str(trialNum), h.targetFlag)
    set(h.totalTrialNumUI,'String',num2str(trialNum));
   
    if trialNum > 1 
        mLatency(trialNum,6) = toc(latencyTic);
    end

    % Trial phases with precision timestamps
    latencyTic = tic;
    trialCue;      
    mLatency(trialNum,1) = toc(latencyTic);

    h.postCuePeriodlatencyTic = tic;
    postCuePeriod; 
    mLatency(trialNum,2) = toc(h.postCuePeriodlatencyTic);

    latencyTic = tic;
    visiStim;      
    mLatency(trialNum,3) = toc(latencyTic);

    latencyTic = tic;
    responseWindow;
    mLatency(trialNum,4) = toc(latencyTic);

    latencyTic =tic;
    ITIperiod;      
    mLatency(trialNum,5) = toc(latencyTic);
    
    latencyTic =tic;

    % 计算lick rate数据（保留数据收集但不绘制热图）
    calculateLickRateData(trialNum);

    fprintf('In this trial, lick in RW = %s times! ',num2str(h.lickInRWOneTrial));
    fprintf('lick out of RW = %s times! \n',num2str(h.outRWCounterSingle));
    
    h.correctRate = (h.hitTrialNumber + h.CRTrialNumber) / trialNum;
    
    % 扩展数据记录
    h.data1(trialNum,1) = trialNum;
    h.data1(trialNum,2) = h.resultFlag;   
    h.data1(trialNum,3) = h.lickInRWOneTrial;
    h.data1(trialNum,4) = h.outRWCounterSingle;
    h.data1(trialNum,5) = round(toc(h.trialGlobalTic),4); % 提高精度到0.1ms
    h.data1(trialNum,6) = round(h.correctRate,4);
    h.data1(trialNum,7) = h.targetFlag; % 记录target/non-target
    h.data1(trialNum,8) = length(h.trialLickTimes); % 记录总lick次数
    
    plot(h.ratePlot,trialNum,h.correctRate,'-ok');
    
    % 重置指示器
    if h.targetFlag     
       set(h.targetBox,'Visible','off');
     elseif ~h.targetFlag
        set(h.nontargetBox,'Visible','off');
    end

    % 更新UI背景色
    switch h.resultFlag
        case 1
         set(h.hitCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
        case 2
        set(h.missCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
        case 3
        set(h.FACounterUI,'BackgroundColor',[0.9 0.9 0.9]);
        case 4
        set(h.CRCounterUI,'BackgroundColor',[0.9 0.9 0.9]);
    end

end

% 清理和保存 - 包括新增的水泵timer
stop(h.tLickCounter);
delete(h.tLickCounter);
delete(h.tRefractory);
delete(h.tAirpuff);
if exist('h.tWaterPump','var')
    if strcmp(h.tWaterPump.Running,'on')
        stop(h.tWaterPump);
    end
    delete(h.tWaterPump);
end

disp('---------------------finished!---------------')
totalTime = toc(totalTic);
totalLickTimes = 0;
for i = 1 : numel(h.lickdata)
    totalLickTimes = totalLickTimes + length(h.lickdata{i});
end
PsychPortAudio('Close',h.audioHandle); 
fprintf('>> total time cost:  %s minutes %s seconds \n',num2str(floor((totalTime)/60)),num2str(mod(totalTime,60)));
fprintf('>> total lick times: %s times \n',num2str(h.inRWCounter));

% 保存完整数据，包括高精度时间戳
savestr = [datestr(now, 'yyyymmdd') '_'  h.mouseID{1} 'discri.mat'];
save(savestr, 'h', 'mLatency','f');
fprintf('data saved as, %s', savestr);
sca

%% functions - 整合stage2和stage3的所有函数

% 只保留数据计算，删除热图绘制
function calculateLickRateData(trialNum)
global h
    % 计算lick rate数据但不绘制热图
    lickCount = zeros(size(h.timeAxis));
    
    for i = 1:length(h.timeAxis)
        windowStart = h.timeAxis(i) - h.lickRateTimeWindow/2;
        windowEnd = h.timeAxis(i) + h.lickRateTimeWindow/2;
        
        licksInWindow = sum(h.trialLickTimes >= windowStart & h.trialLickTimes <= windowEnd);
        lickCount(i) = licksInWindow; 
    end
    
    % 保存数据供后续分析使用
    h.lickRateTimeCourse{trialNum} = lickCount;
    h.lickRateMatrix(trialNum, :) = lickCount;
end

function infoUI(~,~)
 global h
 prompt = {'mouseID', 'trainStage', 'dayNumber', 'saveDir'};
 dlgtitle = 'mouse information';
 dims = [1 35];
 definput = {'', 'stage3', '', '/home/liu/github/original_program/behavioral_program/data'};
 h.mouseID = inputdlg(prompt, dlgtitle, dims, definput);
end

function updateVbl(~,~)
 global h
 h.vbl = Screen('Flip',h.window);
 h.vblt0 = h.vbl;
end

function countDown(~,~)
 countdownDuration = 180; % 保持stage3的3分钟倒计时
 disp('---------Countdown started...check the mouse and lick spout!!!--------------------')

 for remainingSeconds = countdownDuration :-1 :0
     fprintf('Time remaining: %d seconds \n',remainingSeconds);
     pause(1);
 end
 disp('Countdown finished! Start conditioning!')
end

%% trial procedure functions - 保留stage3功能但添加高精度时间戳
function trialCue(~,~)
global h
PsychPortAudio('Start',h.audioHandle,1,0,1);
[startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop',h.audioHandle,1,1);
disp('>>Trial Cue ended!! Post-cue period starts!!')
end

function postCuePeriod(~,~)
global h
    h.inOrOutRW = -1;  
    
    if ~strcmp(h.tLickCounter.Running,'on')
        start(h.tLickCounter); 
    end
    
    h.postCueTime = tic;
    while toc(h.postCueTime) <= 1 
        WaitSecs(0.001); % 减少CPU占用
    end
    disp('>>post-cue period finished! Visual stimulation starts!')
end

function visiStim(~,~)
  global h
  h.inOrOutRW = 0; 
  
  % 使用stage2的优化视觉刺激代码
  Screen('FillRect', h.window, h.grey);
  vbl = Screen('Flip', h.window);
  vblt0 = vbl;
  h.vbl = vbl; h.vblt0 = vblt0; % 为兼容性保留
  startTime = vbl;
  
  h.propertiesMat(1) = 0; % 重置相位
  
  VSLength = 1; 
  frameCount = 0;
  
  while (vbl - vblt0) <= VSLength
    frameCount = frameCount + 1;
    
    Screen('DrawTextures', h.window, h.gabortex, [], [], h.visiOri, [], [], [], [],...
          kPsychDontDoRotation, h.propertiesMat');
    
    nextFlipTime = startTime + frameCount * h.ifi;
    vbl = Screen('Flip', h.window, nextFlipTime - 0.5 * h.ifi);
    h.vbl = vbl; % 更新全局vbl以兼容lick检测
    
    h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
  end
  
  Screen('FillRect', h.window, h.grey);
  Screen('Flip', h.window);
  
  h.inOrOutRW = 1;
  disp('>>Visual stimulation ended!! Response window starts!!')
end

function responseWindow(~,~)
global h
 h.rwTime = tic;
 h.rwLimit = 4;
  
  while toc(h.rwTime) <= h.rwLimit    
      WaitSecs(0.001);
  end
  
  if h.lickInRWOneTrial > 0
      h.licktrial = h.licktrial + 1;
  end
  
h.inOrOutRW = 0;   
disp('>>Response window ended!! ITI start!')
end

function ITIperiod(~,~)
global h
    h.inOrOutRW = 2;
    
    % 判断trial结果 - 保留stage3的完整逻辑
    if h.lickInRWOneTrial == 0   
        if h.targetFlag    
            h.missTrialNumber = h.missTrialNumber + 1;
            h.resultFlag = 2; % Miss
            set(h.missTrialNumUI,'String',num2str(h.missTrialNumber));
            set(h.missCounterUI,'BackgroundColor',[1 0 0]);
        elseif ~h.targetFlag   
            h.CRTrialNumber = h.CRTrialNumber + 1;
            h.resultFlag = 4; % Correct Rejection
            set(h.CRTrialNumUI,'String',num2str(h.CRTrialNumber));
            set(h.CRCounterUI,'BackgroundColor',[0 1 0]);
        end
    end
    
    h.ITITime = tic;
    itiDuration = h.ITIperiod(h.trialNum);
    while toc(h.ITITime) <= itiDuration 
        WaitSecs(0.01);
    end
    
    fprintf('>> ITI period finished!!! Time length is %s \n',num2str(itiDuration))
end

function pinStatusChanged(~,~)
global h
    % 安全检查
    if ~strcmp(h.tLickCounter.Running,'on')
        return;
    end
    
    RWflag = h.inOrOutRW;
    lickFlag = true;
    trialFlag = true;
    lickTimesReporter = 0;
    pinValue = readDigitalPin(h.a,h.sensorPin);
    
    if pinValue == true && lickFlag == true
        % 记录高精度lick时间戳
        trialElapsedTime = toc(h.trialGlobalTic);
        h.trialLickTimes = [h.trialLickTimes, trialElapsedTime];
        
        switch RWflag
            case 1 % RW - 保留stage3的完整奖励/惩罚逻辑
                h.inRWCounter = h.inRWCounter + 1;
              if h.targetFlag  % Target trial
                if h.lickInRWOneTrial < 1      % Hit
                    if round(toc(h.rwTime),3) > 10
                        plot(h.trialRaster,1.0005,h.trialNum,'.g');  
                    else
                        plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.g'); 
                    end
                   % 改进的水泵控制 - 确保可靠性
                   try
                       % 先确保水泵处于关闭状态
                       writeDigitalPin(h.a,'D9',0);
                       pause(0.01); % 短暂延时确保状态清除
                       % 开启水泵
                       writeDigitalPin(h.a,'D9',1);
                       % 启动专用timer控制水泵持续时间
                       if strcmp(h.tWaterPump.Running,'off')
                           start(h.tWaterPump);
                       end
                       fprintf('Water pump activated for Hit! ');
                   catch ME
                       warning('Water pump control error: %s', ME.message);
                       % 备用控制方案
                       writeDigitalPin(h.a,'D9',1);
                       pause(0.15);
                       writeDigitalPin(h.a,'D9',0);
                   end
                   h.hitTrialNumber = h.hitTrialNumber + 1; 
                   h.resultFlag = 1;
                   set(h.hitTrialNumUI,'String',num2str(h.hitTrialNumber));
                   set(h.hitCounterUI,'BackgroundColor',[0 1 0]);
                else
                    plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.','Color',[0.6 0.6 0.6]); 
                end
                fprintf('lick in RW %s ',num2str(h.inRWCounter));
                start(h.tRefractory);
                lickFlag = false; 
                
              elseif ~h.targetFlag  % Non-target trial
                if h.lickInRWOneTrial < 1     % False Alarm
                     if round(toc(h.rwTime),3) > 10
                        plot(h.trialRaster,1.001,h.trialNum,'.r');  
                    else
                        plot(h.trialRaster,round(toc(h.rwTime),3) + 1, h.trialNum,'.r');
                    end
                    writeDigitalPin(h.a,'D3',1); % Air puff
                    h.FATrialNumber = h.FATrialNumber + 1;
                    h.resultFlag = 3;
                    set(h.FATrialNumUI,'String',num2str(h.FATrialNumber));   
                    set(h.FACounterUI,'BackgroundColor',[1 0 0]);                                    
                    start(h.tAirpuff);
                    lickFlag = false; 
                    disp('FA licking! Air puff and time-out starts!')
               else
                    plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.','Color',[0.6 0.6 0.6]); 
                end
                    start(h.tRefractory);
              end

              % 记录lick数据
              if lickFlag == false
                if trialFlag == true
                   trialFlag = false;
                   fprintf('licktrial = %s ',num2str(h.licktrial));
                   set(h.lickTrialNumUI,'String',num2str(h.licktrial));
                end
                lickTimesReporter = lickTimesReporter + 1;
                h.lickInRWOneTrial = h.lickInRWOneTrial + lickTimesReporter;
                tocReporter = toc(h.rwTime);
                fprintf('licktime = %s  @%s seconds \n',num2str(h.lickInRWOneTrial),num2str(round(tocReporter,3)));
                h.lickdata{h.licktrial}(h.lickInRWOneTrial,1) = h.lickInRWOneTrial;
                h.lickdata{h.licktrial}(h.lickInRWOneTrial,2) = tocReporter;
              end

            case -1   % post-cue period                                                                        
                    h.postCueTime = tic;
                    disp('!!Early lick! Reset post cue period timer!!')
                    start(h.tRefractory);
            case 0    % visual stimulation period
                 if exist('h.vbl','var') && exist('h.vblt0','var') && (h.vbl - h.vblt0) <= 1
                    plot(h.trialRaster,round((h.vbl - h.vblt0),3),h.trialNum,'.','Color',[0.6 0.6 0.6]); 
                 end 

            case 2   % ITI
                plot(h.trialRaster,round(toc(h.ITITime),3)+1+h.rwLimit,h.trialNum,'.','Color',[0.6 0.6 0.6]); 
          h.outRWCounter = h.outRWCounter + 1;
          h.outRWCounterSingle = h.outRWCounterSingle + 1;
          start(h.tRefractory);
        end
   end
end

% Timer callback functions - 保留stage3的所有timer功能
function refStart(~,~)
    global h
      if strcmp(h.tLickCounter.Running,'on')
          stop(h.tLickCounter);
      end
      % 确保水泵关闭 - 改进的控制逻辑
      try
          writeDigitalPin(h.a,'D9',0);
          % 停止水泵timer如果正在运行
          if strcmp(h.tWaterPump.Running,'on')
              stop(h.tWaterPump);
          end
      catch ME
          warning('Water pump shutdown error in refStart: %s', ME.message);
      end
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
 writeDigitalPin(h.a,'D3',0);
 h.rwLimit = toc(h.rwTime) + 7; % 延长RW时间作为惩罚
 stop(h.tAirpuff);
end

function waterPumpEnd(~,~)
 global h
 % 水泵专用关闭函数 - 确保可靠关闭
 try
     writeDigitalPin(h.a,'D9',0);
     fprintf('Water pump deactivated. ');
 catch ME
     warning('Water pump deactivation error: %s', ME.message);
     % 尝试多次关闭
     for i = 1:3
         try
             pause(0.01);
             writeDigitalPin(h.a,'D9',0);
             break;
         catch
             continue;
         end
     end
 end
 
 if strcmp(h.tWaterPump.Running,'on')
     stop(h.tWaterPump);
 end
end

function seq = randomSequence(n, m)
seq = zeros(n, 1);
seq(1 : m) = 1;
seq(randperm(n)) = seq;
nPermit = 3;
iSame = 0;
  for i = 1 : n - 1
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
