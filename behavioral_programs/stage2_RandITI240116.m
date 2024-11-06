%% modified program from stageTwoRasterPlot.m 24.1.5
% new functions
%1. randomized post-cue period following normal distribution 
% with mean=1, SD = 0.1
%2. modified ploting tic/toc timer setting
timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor, h the gray color background
global h 
PsychDefaultSetup(2);
%Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);
h.screenNumber = max(Screen('Screens'));
h.white = WhiteIndex(h.screenNumber);
h.grey = h.white / 2;
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution); 
h.ifi = Screen('GetFlipInterval',h.window); 
h.topPriorityLevel = MaxPriority(h.window); 
Priority(h.topPriorityLevel);

%% initialize sound configuration
InitializePsychSound;

%open psych-audio port
h.sampleF = 48000;
h.audioHandle = PsychPortAudio('Open',12, 1, 1, h.sampleF, 2);  %use PsychPortAudio('GetDevices') to find steinberg UR12, and change the first number with UR12 number
PsychPortAudio('Volume', h.audioHandle, 0.02);      %auditory cue volume

%pre-allocate audio buffer
[myBeep, samplingRate] = MakeBeep(10000, 0.1, h.sampleF);
buffer = [myBeep;myBeep];
PsychPortAudio('FillBuffer',h.audioHandle,buffer);

%% start communication
h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
h.sensorPin = 'D13';
h.waterPumpPin = 'D9';

%% Gabor presetting
%size of the gabor patch, full of the height in this case
h.gaborDimPix = h.windowRect(4)*2;
h.width = h.windowRect(3);
h.height = h.windowRect(4);

%center of diplay position
h.center = [(h.width-h.height)/2,0,(h.width+h.height)/2,h.height];

%other parameters
h.sigma = h.gaborDimPix;
h.orientationTarget = 0; %vertical
h.orientationNontarget = 90; %horizontal
h.contrast = 1;
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
h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];
updateVbl;
h.waitframes = 1;
h.phasePerFrame = 5 * pi;  %change the speed of grating moving

%% program variables
trialNumLimit = 300;
h.licktrial = 1;
lickTrialLimit = 150;
h.lickdata = {};
h.data1 = zeros(trialNumLimit,4);
h.totalLickTimes = 0;

%create random ITI time length ranging from 4 ~ 6s.
h.ITIperiod = 4 + rand([1 trialNumLimit])*2;
%% display UI, waiting for the initialization
infoUI();

%% counterUI
screenSize = get(0,'Screensize'); screenSize(3) = screenSize(3)/2;       
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);          %open the lick monitor UI

%total counter UI
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.05 0.95 0.1 0.04],'String','total trial number','BackgroundColor',[0 1 1],...
    'FontSize',14);
 h.totalTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.05 0.91 0.1 0.03],'FontSize',14); 
 uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.18 0.95 0.1 0.04],'String','lick trial number','BackgroundColor',[0 1 1],...
    'FontSize',14);
 h.lickTrialNumUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.18 0.91 0.1 0.03],'FontSize',14); 
 % function of hit trial number
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.3 0.95 0.1 0.04],'String','Hit trial number','BackgroundColor',[1 1 0],...
    'FontSize',14);
 h.hitTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.3 0.91 0.1 0.03],'FontSize',14); 

 %visual stimulation parameters
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.11 0.86 0.25 0.03],'String','duration info','BackgroundColor',[0.83 0.5 1],...
     'FontSize',14);

 VSPara = {'VS' 'RW ' 'ITI'};
 VSParaValue = [1 4 4];
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.05 0.81 0.1 0.04],'String',VSPara{1},'BackgroundColor',[1 1 1],...
     'FontSize',14);
 h.tempF =uicontrol(f,'Style','edit','String',num2str(VSParaValue(1)),'Units','normalized',...
    'Position',[0.05 0.79 0.11 0.03],'FontSize',14); 
 
 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.18 0.81 0.1 0.04],'String',VSPara{2},'BackgroundColor',[1 1 1],...
     'FontSize',14);
 h.spatF =uicontrol(f,'Style','edit','String',num2str(VSParaValue(2)),'Units','normalized',...
    'Position',[0.18 0.79 0.11 0.03],'FontSize',14); 

 uicontrol(f,'Style','text','Units','normalized',...
     'Position',[0.3 0.81 0.1 0.04],'String',VSPara{3},'BackgroundColor',[1 1 1],...
     'FontSize',14);
 h.duration =uicontrol(f,'Style','edit','String',num2str(VSParaValue(3)),'Units','normalized',...
    'Position',[0.3 0.79 0.11 0.03],'FontSize',14); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.5 0.95 0.15 0.04],'String','lick in RW times   TOTAL','BackgroundColor',[1 1 1],...
    'FontSize',14);
h.lickInRWtotal = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.5 0.9 0.15 0.04],'FontSize',14); %replace the 0 with callback function

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.7 0.95 0.15 0.04],'String','lick out RW times  TOTAL','BackgroundColor',[1 1 1],...
    'FontSize',14);
 h.lickOutRWtotal = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.7 0.9 0.15 0.04],'FontSize',14); %replace the 0 with callback function

%single trial counter UI
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.5 0.83 0.15 0.04],'String','lick in RW times   SINGLE TRIAL','BackgroundColor',[0 1 0],...
    'FontSize',14);
h.lickInRWsingle = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.5 0.78 0.15 0.04],'FontSize',14); %replace the 0 with callback function

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.7 0.83 0.15 0.04],'String','lick out RW times  SINGLE TRIAL','BackgroundColor',[1 0 0],...
    'FontSize',14);
h.lickOutRWSingle = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.7 0.78 0.15 0.04],'FontSize',14); %replace the 0 with callback function

% trial raster
h.trialRaster = axes(f,'Position',[0.1 0.44 0.3 0.3],'FontSize',14);
title('trial raster');
xlabel('seconds');
ylabel('trial number');
xlim([0 11]);
ylim([0 300]);
set(h.trialRaster,'Ydir','reverse');
hold on

% time length monitor
h.trialTime = axes(f,'Position',[0.55 0.44 0.3 0.3],'FontSize',14);
title('trial time length');
xlabel('trial number');     
xlim([1 300]);
ylabel('seconds');
hold on

%open the back camera
h.backCamUI = axes(f,'Position',[0.55 0.05 0.35 0.35],'Title','Back Camera');
h.backCam = webcam(1);
h.backCamSize = str2double(strsplit(h.backCam.Resolution,'x'));
h.im = image(zeros(h.backCamSize),'Parent',h.backCamUI);
preview(h.backCam,h.im);

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
h.outRWCounter = 0;
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
h.hitTrialNumber = 0;
h.earlyTrialNumber = 0;
h.trialNum = [];
h.rwTime = [];
mLatency = zeros(trialNumLimit,6);   %latency test
%% main loop
start(h.tLickCounter);   %start the lick counter
for trialNum = 1:1000
       h.trialNum = trialNum;
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end                                             %stop the tRefractory
stop(h.tLickCounter)
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
set(h.lickInRWsingle,'String',num2str(h.lickInRWOneTrial));
set(h.lickOutRWSingle,'String',num2str(h.outRWCounterSingle));

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

    fprintf('In this trial, lick in RW = %s times! ',num2str(h.lickInRWOneTrial));
    fprintf('lick out of RW = %s times! \n',num2str(h.outRWCounterSingle));
    h.data1(trialNum,1) = trialNum;
    h.data1(trialNum,2) = h.lickInRWOneTrial;
    h.data1(trialNum,3) = h.outRWCounterSingle;
    h.data1(trialNum,4) = toc(h.trialGlobalTic); %time length for individual trials
    trialTimeLength = h.data1(trialNum,4);
    plot(h.trialTime,h.data1(trialNum,1),h.data1(trialNum,4),'-ok');
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
h.lickdata
h.data1
save h

%% functions
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
 countdownDuration = 180;
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
    start(h.tLickCounter); 
    h.postCueTime = tic;
    while toc(h.postCueTime) <=  1     %fixed 1 seconds
    end
    disp('>>post-cue period finished! Visual stimulation starts!')
end

function visiStim(~,~)
  global h
  h.inOrOutRW = 0; 
  updateVbl;
  VSLength = 1; %visual stimulation in the first second
  while h.vbl - h.vblt0 <= VSLength
  Screen('DrawTextures', h.window, h.gabortex, [], [], h.orientationTarget, [], [], [], [],...
        kPsychDontDoRotation, h.propertiesMat');
  h.vbl = Screen('Flip', h.window, h.vbl + (h.waitframes - 0.5) * h.ifi);   %double buffering 
  h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
%       if h.vbl - h.vblt0 > VSLength*0.9
%           h.inOrOutRW = 1;   %the last 10% of the Visi stim is regarded as RW
%       end
  end
  updateVbl;
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
    while toc(h.ITITime) <= h.ITIperiod(h.trialNum)    %index the randomized ITI period
    end
    fprintf('>> ITI period finished!!! Time length is %s \n',num2str(h.ITIperiod(h.trialNum)))
end

function pinStatusChanged(~,~)
global h
    RWflag = h.inOrOutRW;
    lickFlag = true;
    trialFlag = true;
    lickTimesReporter = 0;
    pinValue = readDigitalPin(h.a,h.sensorPin);
    if pinValue == true && lickFlag == true
        switch RWflag
            case 1 % lick in RW
                h.inRWCounter = h.inRWCounter + 1;
                set(h.lickInRWtotal,'String',num2str(h.inRWCounter));
                if h.lickInRWOneTrial < 1 %&& h.outRWCounterSingle == 0  %only the first  lick is rewarded
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
             set(h.lickInRWsingle,'String',num2str(h.lickInRWOneTrial));
             fprintf('licktime = %s  @%s seconds \n',num2str(h.lickInRWOneTrial),num2str(tocReporter));
             h.lickdata{h.licktrial}(h.lickInRWOneTrial,1) = h.lickInRWOneTrial;
             h.lickdata{h.licktrial}(h.lickInRWOneTrial,2) = tocReporter;
           end
            case -1                                                                            % early lick reset the post-cue period timer
                    h.postCueTime = tic;
                    disp('!!Early lick! Reset post cue period timer!!')
                    start(h.tRefractory);
            otherwise
                   if RWflag == 0 && (h.vbl - h.vblt0) <= 1
                       plot(h.trialRaster,round((h.vbl - h.vblt0),3),h.trialNum,'.','Color',[0.6 0.6 0.6]);  %plot the early lick in visual stimulation window
                   elseif RWflag == 2
                       plot(h.trialRaster,round(toc(h.ITITime),3)+5,h.trialNum,'.','Color',[0.6 0.6 0.6]); %plot the ITI licking in 5~9 s ITI window
                   end
                h.outRWCounter = h.outRWCounter + 1;
                h.outRWCounterSingle = h.outRWCounterSingle + 1;
                set(h.lickOutRWtotal,'String',num2str(h.outRWCounter));
                set(h.lickOutRWSingle,'String',num2str(h.outRWCounterSingle));
                fprintf('Lick out of RW in this trial: %s  total: %s \n',...
                    num2str(h.outRWCounterSingle),num2str(h.outRWCounter));  
                start(h.tRefractory);
        end
    end
end

function refStart(~,~)
    global h
      stop(h.tLickCounter);
      %disp('timer 1 stopped by timer-2 StartFcn');
      writeDigitalPin(h.a,'D9',0);
end

function refEnd(~,~)
    global h
    start(h.tLickCounter);
    stop(h.tRefractory);
end
