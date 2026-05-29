%% modified program from stageThreeDemo240129.m
%major changes at stageFour:
%1: 150 trials in total, 100%, 25%, 6.25% trials are 50 for each contrast.
%2: 2 sessions for imaging everyday
%3: separated index for contrast and target/non-target
timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor, h the gray color background
global h 
PsychDefaultSetup(2);

% 添加性能优化的Screen偏好设置（来自stage3）
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
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution); 
h.ifi = Screen('GetFlipInterval',h.window); 
h.topPriorityLevel = MaxPriority(h.window); 
Priority(h.topPriorityLevel);

% 预热显卡和优化缓存（来自stage3）
for i = 1:5
    Screen('Flip', h.window);
end

%% initialize sound configuration
InitializePsychSound;

%open psych-audio port - 使用stage3的设备查找方式
h.sampleF = 48000;
deviceList = PsychPortAudio('GetDevices');

% find the device ID for the Steinberg UR12（来自stage3）
device_id = [];
for i = 1:length(deviceList)
    if contains(deviceList(i).DeviceName, 'UR12') || contains(deviceList(i).DeviceName, 'Steinberg')
        device_id = i;
        disp(['Audio device found with Device ID: ' num2str(device_id)]);
        break;
    end
end
if isempty(device_id)
    device_id = 13; % 回退方案
    warning('Steinberg UR12 device not found. Using fallback device ID: %d', device_id);
end

h.audioHandle = PsychPortAudio('Open', device_id-1, 1, 1, h.sampleF, 2);
PsychPortAudio('Volume', h.audioHandle, 0.02);      %auditory cue volume 

%pre-allocate audio buffer
[myBeep, samplingRate] = MakeBeep(10000, 0.1, h.sampleF);
buffer = [myBeep;myBeep];
PsychPortAudio('FillBuffer',h.audioHandle,buffer);

%% start communication
h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
h.sensorPin = 'D13';
h.waterPumpPin = 'D9';
h.airPumpPin = 'D3';

%% Gabor presetting
%size of the gabor patch, full of the height in this case
h.gaborDimPix = h.windowRect(4)*2;
h.width = h.windowRect(3);
h.height = h.windowRect(4);

%center of diplay position
h.center = [(h.width-h.height)/2,0,(h.width+h.height)/2,h.height];

%other parameters
h.sigma = h.gaborDimPix;
h.oriTarget = 0; %vertical
h.oriNonTarget = 90; %horizontal
h.contrast = [1 0.1 0.01 0.001];              % 100%, 10%, 1%, 0.1%
h.aspectRatio = 1;
h.phase = 0; 

%spatial frequency
h.numCycles = 7;
h.freq = h.numCycles / h.gaborDimPix;

%make procedural gabor texture
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;
h.gabortex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

%make the property matrix

h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast(1), h.aspectRatio, 0, 0, 0];

updateVbl;
h.waitframes = 1;
h.phasePerFrame = 5 * pi;  %change the speed of grating moving

% 预计算常用值以提高性能（来自stage3）
h.halfIfi = 0.5 * h.ifi;
h.frameAdvance = (h.waitframes - 0.5) * h.ifi;

%% program variables
nTrial = 160;  %total 150 trials in a block, tar vs nontar is 50/50
h.tarAndContrastCombo = [h.oriTarget h.contrast(1);...
                                                   h.oriTarget h.contrast(2);...
                                                   h.oriTarget h.contrast(3);...
                                                   h.oriTarget h.contrast(4);...
                                                   h.oriNonTarget h.contrast(1);...
                                                   h.oriNonTarget h.contrast(2);
                                                   h.oriNonTarget h.contrast(3);...
                                                   h.oriNonTarget h.contrast(4)];   % 6 combination of target/nontarget and contrast
h.randSeq = randomSequence(nTrial,h.tarAndContrastCombo);  %randomize the target/nontarget sequence
h.oriSequence = h.randSeq(:,1);    %the randomized orientation sequence
h.contrastSequence = h.randSeq(:,2);
h.licktrial = 1;
h.lickdata = {};
h.data1 = zeros(nTrial,10); % 扩展为10列：trial#, result, lickInRW, lickOutRW, trialTime, correctRate, contrast, targetFlag, totalLicks, firstLickLatency

%trial counters for UI
h.totalLickTimes = 0;
h.hitTrialNumber = 0;
h.FATrialNumber = 0;
h.CRTrialNumber = 0;
h.missTrialNumber = 0;
h.resultFlag = [];
h.correctRate = 0;

% Lick rate monitoring variables - 引入stage3的功能
h.lickRateTimeWindow = 0.1; 
h.lickRateTimeStep = 0.05; 
h.trialLickTimes = []; 
h.lickRateTimeCourse = {}; 
h.timeAxis = 0:h.lickRateTimeStep:15; % 扩展到15秒

% 热图数据矩阵
h.lickRateMatrix = zeros(nTrial, length(h.timeAxis));
h.heatmapHandle = [];

%create random ITI time length ranging from 4 ~ 6s.
h.ITIperiod = 4 + rand([1 nTrial])*2;

%% UI, waiting for the initialization
infoUI();

%% UI - 使用stage3的优化UI布局
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);          %open the lick monitor UI                        

%total counter UI
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.05 0.95 0.1 0.04],'String','total trial number','BackgroundColor',[0 1 1],...
    'FontSize',16);
 h.totalTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.05 0.91 0.1 0.03],'FontSize',16); 
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.18 0.95 0.1 0.04],'String','lick trial number','BackgroundColor',[0 1 1],...
    'FontSize',16);
 h.lickTrialNumUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.18 0.91 0.1 0.03],'FontSize',16); 

 % function of result counter
h.hitCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.5 0.95 0.15 0.04],'String','Hit trial number','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',16);
 h.hitTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.5 0.91 0.15 0.03],'FontSize',16); 

 h.missCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.7 0.95 0.15 0.04],'String','Miss trial number','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',16);
 h.missTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.7 0.91 0.15 0.03],'FontSize',16); 

 h.FACounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.5 0.83 0.15 0.04],'String','FA trial number','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',16);
 h.FATrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.5 0.79 0.15 0.03],'FontSize',16); 

 h.CRCounterUI = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.7 0.83 0.15 0.04],'String','CR trial number','BackgroundColor',[0.9 0.9 0.9],...
    'FontSize',16);
 h.CRTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.7 0.79 0.15 0.03],'FontSize',16); 

 %target and nontarget box indicator
 h.targetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.93 0.1 0.05],'String','Target trial starts','BackgroundColor',[1 1 0],...
    'FontSize',16,'Visible','off');
 h.nontargetBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.82 0.1 0.05],'String','Nontarget trial starts','BackgroundColor',[1 1 0],...
    'FontSize',16,'Visible','off');

 %contrast box indicator
  h.hundredPerContrastBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.68 0.1 0.05],'String','100% contrast','BackgroundColor',[0 0.7 0.7],...
    'FontSize',16,'Visible','off');

    h.tenPerContrastBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.6 0.1 0.05],'String','10% contrast','BackgroundColor',[0 0.7 0.7],...
    'FontSize',16,'Visible','off');

     h.onePerContrastBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.52 0.1 0.05],'String','1% contrast','BackgroundColor',[0 0.7 0.7],...
    'FontSize',16,'Visible','off');

     h.oneTenthContrastBox = uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.88 0.44 0.1 0.05],'String','0.1% contrast','BackgroundColor',[0 0.7 0.7],...
    'FontSize',16,'Visible','off');
     
 %visual stimulation parameters
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.11 0.86 0.25 0.03],'String','duration info','BackgroundColor',[0.83 0.5 1],...
     'FontSize',16);

 VSPara = {'VS' 'RW ' 'ITI'};
 VSParaValue = [1 4 4];
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.05 0.81 0.1 0.04],'String',VSPara{1},'BackgroundColor',[1 1 1],...
     'FontSize',16);
 h.tempF =uicontrol(f,'Style','edit','String',num2str(VSParaValue(1)),'Units','normalized',...
    'Position',[0.05 0.79 0.11 0.03],'FontSize',16); 
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.18 0.81 0.1 0.04],'String',VSPara{2},'BackgroundColor',[1 1 1],...
     'FontSize',16);
 h.spatF =uicontrol(f,'Style','edit','String',num2str(VSParaValue(2)),'Units','normalized',...
    'Position',[0.18 0.79 0.11 0.03],'FontSize',16); 

 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.3 0.81 0.1 0.04],'String',VSPara{3},'BackgroundColor',[1 1 1],...
     'FontSize',16);
 h.duration =uicontrol(f,'Style','edit','String','4~6','Units','normalized',...
    'Position',[0.3 0.79 0.11 0.03],'FontSize',16); 

% trial raster
h.trialRaster = axes(f,'Position',[0.1 0.44 0.3 0.3]);
title('trial raster');
xlabel('seconds');
ylabel('trial number');
xlim([0 15]);
ylim([1 nTrial]);
set(h.trialRaster,'Ydir','reverse');
hold on

% success rate monitor
h.ratePlot = axes(f,'Position',[0.55 0.44 0.3 0.3]);
title('correct rate plot');
xlabel('trial number');     
xlim([1 nTrial]);
ylim([0 1]);
ylabel('correct rate');
hold on

%% open cameras - 使用stage3的双摄像头功能
try
    h.backCam = webcam(1);
    h.frontCam = webcam(2);

    backRes = str2double(strsplit(h.backCam.Resolution,'x'));
    frontRes = str2double(strsplit(h.frontCam.Resolution,'x'));

    % Calculate aspect ratios
    backAspect = backRes(1) / backRes(2);
    frontAspect = frontRes(1) / frontRes(2);

    % Set UI positions based on aspect ratios
    backHeight = 0.28;
    backWidth = backAspect * backHeight;
    frontHeight = 0.28;
    frontWidth = frontAspect * frontHeight;

    % Place back camera UI
    h.backCamUI = axes(f, 'Position', [0.55 0.05 backWidth backHeight]);
    h.im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', h.backCamUI);
    preview(h.backCam, h.im);
    text(h.backCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);

    % Place front camera UI
    h.frontCamUI = axes(f, 'Position', [0.1 0.05 frontWidth frontHeight]);
    h.im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', h.frontCamUI);
    preview(h.frontCam, h.im2);
    text(h.frontCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
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

% 新增：水泵控制专用timer，确保奖励持续时间（来自stage3）
h.tWaterPump = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.15,...
       'TimerFcn',@(src,event)waterPumpEnd);

%% 3 minutes countdown 
countDown;

%% start tic 
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
%start(h.tLickCounter);   %start the lick counter
for trialNum = 1:nTrial
       h.trialNum = trialNum;
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end                                      %stop the tRefractory at the beginning of a block
stop(h.tLickCounter);
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
h.trialLickTimes = []; % 重置当前trial的lick时间记录
h.visiOri = h.oriSequence(h.trialNum); 
h.contrastOfThisTrial = h.contrastSequence(h.trialNum);%index the visi orientation using trial number
if h.contrastOfThisTrial == 1     %indicator box of contrast trials
   set(h.hundredPerContrastBox,'Visible','on');
elseif h.contrastOfThisTrial == 0.1
    set(h.tenPerContrastBox,'Visible','on');
elseif h.contrastOfThisTrial == 0.01
     set(h.onePerContrastBox,'Visible','on');
elseif h.contrastOfThisTrial == 0.001
    set(h.oneTenthContrastBox,'Visible','on');
end

h.targetFlag = (h.visiOri == h.oriTarget); % if target, targetFlag is 1; if not target Flag is 0
if h.targetFlag     %indicator box of target/nontarget trials
   set(h.targetBox,'Visible','on');
elseif ~h.targetFlag
    set(h.nontargetBox,'Visible','on');
end

% 记录trial开始的高精度时间戳（来自stage3）
h.trialGlobalTic = tic;
disp('------------------new trial-----------------------')
fprintf('trialNum = %s, Target = %d, Contrast = %.3f\n', num2str(trialNum), h.targetFlag, h.contrastOfThisTrial)
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

  % 计算lick rate数据（来自stage3）
  calculateLickRateData(trialNum);

  fprintf('In this trial, lick in RW = %s times! ',num2str(h.lickInRWOneTrial));
  fprintf('lick out of RW = %s times! \n',num2str(h.outRWCounterSingle));
  h.correctRate = (h.hitTrialNumber + h.CRTrialNumber) / trialNum;
  
  % 扩展数据记录（来自stage3）
  h.data1(trialNum,1) = trialNum;
  h.data1(trialNum,2) = h.resultFlag;   %save the result in the data1 matrix: 1.hit 2.miss 3.FA 4.CR
  h.data1(trialNum,3) = h.lickInRWOneTrial;
  h.data1(trialNum,4) = h.outRWCounterSingle;
  h.data1(trialNum,5) = round(toc(h.trialGlobalTic),4); %time length for individual trials (高精度)
  h.data1(trialNum,6) = round(h.correctRate,4);
  h.data1(trialNum,7) = h.contrastOfThisTrial; %contrast of this trial
  h.data1(trialNum,8) = h.targetFlag; % 记录target/non-target
  h.data1(trialNum,9) = length(h.trialLickTimes); % 记录总lick次数
  % 计算第一次lick的latency（如果有lick的话）
  if ~isempty(h.trialLickTimes)
      h.data1(trialNum,10) = min(h.trialLickTimes); % 第一次lick相对trial开始的时间
  else
      h.data1(trialNum,10) = NaN;
  end
  
  plot(h.ratePlot,trialNum,h.correctRate,'-ok');
  
  if h.targetFlag     %reset target/non-target indicator
   set(h.targetBox,'Visible','off');
 elseif ~h.targetFlag
    set(h.nontargetBox,'Visible','off');
  end

 if h.contrastOfThisTrial == 1     %reset contrast indicator box
   set(h.hundredPerContrastBox,'Visible','off');
elseif h.contrastOfThisTrial == 0.1
    set(h.tenPerContrastBox,'Visible','off');
elseif h.contrastOfThisTrial == 0.01
     set(h.onePerContrastBox,'Visible','off');
elseif h.contrastOfThisTrial == 0.001
    set(h.oneTenthContrastBox,'Visible','off');
end

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
PsychPortAudio('Close',h.audioHandle); %close the audio stimulation port
fprintf('>> total time cost:  %s minutes %s seconds \n',num2str(floor((totalTime)/60)),num2str(mod(totalTime,60)));
fprintf('>> total lick times: %s times \n',num2str(h.inRWCounter));

% 保存完整数据，包括高精度时间戳（来自stage3）
savestr = [datestr(now, 'yyyymmdd') '_'  h.mouseID{1} '_stage4_discri.mat'];
save(savestr, 'h', 'mLatency','f');
fprintf('data saved as, %s\n', savestr);
sca

%% functions
function infoUI(~,~)
 global h
 prompt = {'mouseID', 'trainStage', 'dayNumber', 'saveDir'};
 dlgtitle = 'mouse information';
 dims = [1 35];
 definput = {'', 'stage4', '', '/home/liu/github/original_program/behavioral_program/data'};
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
    if strcmp(h.tLickCounter.Running, 'off')
        start(h.tLickCounter);
    else
        stop(h.tLickCounter); % Stop the timer if it's running
        start(h.tLickCounter); % Start it again
    end
    h.postCueTime = tic;
    while toc(h.postCueTime) <= 1 
        WaitSecs(0.001); % 减少CPU占用（来自stage3）
    end
    disp('>>post-cue period finished! Visual stimulation starts!')
end

function visiStim(~,~)
  global h
  h.inOrOutRW = 0; 
  
  % 使用stage3的优化视觉刺激代码
  h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrastOfThisTrial, h.aspectRatio, 0, 0, 0];
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
 %h.inOrOutRW = 1;  %set as true when RW starts  
  while toc(h.rwTime) <= h.rwLimit    %the RW last for 4 seconds
      WaitSecs(0.001); % 减少CPU占用（来自stage3）
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
    if h.lickInRWOneTrial == 0   %no lick in RW
        if h.targetFlag    %Miss the target stimulation-->Miss
            h.missTrialNumber = h.missTrialNumber + 1;
            h.resultFlag = 2;
            set(h.missTrialNumUI,'String',num2str(h.missTrialNumber));
             set(h.missCounterUI,'BackgroundColor',[1 0 0]);
        elseif ~h.targetFlag   %correctly reject the non-target stimulation
            h.CRTrialNumber = h.CRTrialNumber + 1;
            h.resultFlag = 4;
            set(h.CRTrialNumUI,'String',num2str(h.CRTrialNumber));
             set(h.CRCounterUI,'BackgroundColor',[0 1 0]);
        end
    end
    h.ITITime = tic;
    itiDuration = h.ITIperiod(h.trialNum);
    while toc(h.ITITime) <= itiDuration %index the randomized ITI period
        WaitSecs(0.01); % 减少CPU占用（来自stage3）
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
        % 记录高精度lick时间戳（来自stage3）
        trialElapsedTime = toc(h.trialGlobalTic);
        h.trialLickTimes = [h.trialLickTimes, trialElapsedTime];
        
        switch RWflag
            case 1 % RW
                h.inRWCounter = h.inRWCounter + 1;
                %set(h.lickInRWtotal,'String',num2str(h.inRWCounter));
              if h.targetFlag  %when target stimulation
                if h.lickInRWOneTrial < 1      %----------------------hit------------------------
                    if round(toc(h.rwTime),3) > 10
                        plot(h.trialRaster,1.0005,h.trialNum,'.g');  
                    else
                        plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.g'); 
                    end
                   % 改进的水泵控制 - 确保可靠性（来自stage3）
                   try
                       writeDigitalPin(h.a,'D9',0);
                       pause(0.01);
                       writeDigitalPin(h.a,'D9',1);
                       if strcmp(h.tWaterPump.Running,'off')
                           start(h.tWaterPump);
                       end
                       fprintf('Water pump activated for Hit! ');
                   catch ME
                       warning('Water pump control error: %s', ME.message);
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
                %trigger the counting and recording system below
              elseif ~h.targetFlag  %when non-target stimulation
                if h.lickInRWOneTrial < 1     %----------------------FA-------------------------
                     if round(toc(h.rwTime),3) > 10
                        plot(h.trialRaster,1.001,h.trialNum,'.r');  
                    else
                        plot(h.trialRaster,round(toc(h.rwTime),3) + 1, h.trialNum,'.r');
                    end
                    writeDigitalPin(h.a,'D3',1);
                    h.FATrialNumber = h.FATrialNumber + 1;
                    h.resultFlag = 3;
                    set(h.FATrialNumUI,'String',num2str(h.FATrialNumber));   
                    set(h.FACounterUI,'BackgroundColor',[1 0 0]);                                    %FA counter UI
                    start(h.tAirpuff);
                    lickFlag = false; 
                    disp('FA licking! Air puff and time-out starts!')
               else
                    plot(h.trialRaster,round(toc(h.rwTime),3) + 1,h.trialNum,'.','Color',[0.6 0.6 0.6]); 
                end
                    %fprintf('lick in RW %s ',num2str(h.inRWCounter));
                    start(h.tRefractory);
              end

              %parallel structure, working when the lickFlag is changed.
              if lickFlag == false
                if trialFlag == true
                   trialFlag = false;
                   fprintf('licktrial = %s ',num2str(h.licktrial));
                   set(h.lickTrialNumUI,'String',num2str(h.licktrial));
                end
                lickTimesReporter = lickTimesReporter + 1;
                h.lickInRWOneTrial = h.lickInRWOneTrial + lickTimesReporter;
                tocReporter = toc(h.rwTime);
                %set(h.lickInRWsingle,'String',num2str(h.lickInRWOneTrial));
                fprintf('licktime = %s  @%s seconds \n',num2str(h.lickInRWOneTrial),num2str(round(tocReporter,3)));
                h.lickdata{h.licktrial}(h.lickInRWOneTrial,1) = h.lickInRWOneTrial;
                h.lickdata{h.licktrial}(h.lickInRWOneTrial,2) = tocReporter;
              end

            case -1   % post-cue period, reset if licked                                                                        
                    h.postCueTime = tic;
                    disp('!!Early lick! Reset post cue period timer!!')
                    start(h.tRefractory);
            case 0    % visual stimulation period
                 if (h.vbl - h.vblt0) <= 1
                    plot(h.trialRaster,round((h.vbl - h.vblt0),3),h.trialNum,'.','Color',[0.6 0.6 0.6]); %plot the early lick in visual stimulation window
                 end 

            case 2   % ITI
                plot(h.trialRaster,round(toc(h.ITITime),3)+1+h.rwLimit,h.trialNum,'.','Color',[0.6 0.6 0.6]); %plot the ITI licking in 5~9 s ITI window
          h.outRWCounter = h.outRWCounter + 1;
          h.outRWCounterSingle = h.outRWCounterSingle + 1;
          %set(h.lickOutRWtotal,'String',num2str(h.outRWCounter));
          %set(h.lickOutRWSingle,'String',num2str(h.outRWCounterSingle));
          %fprintf('Lick out of RW in this trial: %s  total: %s \n',...
             %num2str(h.outRWCounterSingle),num2str(h.outRWCounter));  
          start(h.tRefractory);
        end
   end
end

function refStart(~,~)
    global h
      if strcmp(h.tLickCounter.Running,'on')
          stop(h.tLickCounter);
      end
      % 确保水泵关闭 - 改进的控制逻辑（来自stage3）
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

function airpuffStart(~,~)
  global h
  stop(h.tLickCounter);
  start(h.tRefractory);
end

function airpuffEnd(~,~)
 global h
 writeDigitalPin(h.a,'D3',0);
  h.rwLimit = toc(h.rwTime) + 7;
  stop(h.tAirpuff);
end

function waterPumpEnd(~,~)
 global h
 % 水泵专用关闭函数 - 确保可靠关闭（来自stage3）
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

% 计算lick rate数据（来自stage3）
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

function seq = randomSequence(n, m)
seq = repmat(m,n/length(m),1);
n = size(seq,1);
r = randperm(n);
seq = seq(r,:);
% Up to three consecutive repetitions of the same orientation or contrast are allowed.
nPermit = 3;
iSameA = 0;
iSameB = 0;
  for i = 1 : n - 1
    if seq(i,1) == seq(i + 1,1) %first column
        iSameA = iSameA + 1;
    else
        iSameA = 0;
    end
  
   if seq(i,2) == seq(i + 1,2) %second column
       iSameB = iSameB + 1;
   else
       iSameB = 0;
   end

    j = 0;
    while iSameA > nPermit - 1 || iSameB > nPermit - 1
        temp = seq(end - j,:); 
        seq(end - j,:) = seq(i + 1,:); 
        seq(i + 1,:) = temp; %swap the 4th repetitive element with the end of the sequence

        if seq(i) == seq(i + 1)
            j = j + 1;
        else
            j = 0;
            iSameA = 0;
            iSameB = 0;
        end
    end
  end
end
