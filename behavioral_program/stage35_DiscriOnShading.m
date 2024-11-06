%% modified program from stageTwo.m
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
%Screen('Preference','SkipSyncTests',1);
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

%% initialize sound configuration
InitializePsychSound;

%open psych-audio port
h.sampleF = 48000;
h.audioHandle = PsychPortAudio('Open', 12, 1, 1, h.sampleF, 2);  %use PsychPortAudio('GetDevices') to find UR12, and change the first number with UR12 number
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
nTrial = 160;
nTarget = 80;   %total 300 trials in a block, tar vs nontar is 50/50
h.randSeq = randomSequence(nTrial,nTarget);  %randomize the target/nontarget sequence
h.oriSequence = h.randSeq * 90;    %create the randomized orientation sequence
h.licktrial = 1;
h.lickdata = {};
h.data1 = zeros(nTrial,6);

%trial counters for UI
h.totalLickTimes = 0;
h.hitTrialNumber = 0;
h.FATrialNumber = 0;
h.CRTrialNumber = 0;
h.missTrialNumber = 0;
h.resultFlag = [];
h.correctRate = 0;

%create random ITI time length ranging from 4 ~ 6s.
h.ITIperiod = 4 + rand([1 nTrial])*2;

%% UI, waiting for the initialization
infoUI();

%% UI
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         %
f = figure('Name','trial monitor','Position',screenSize);          %open the lick monitor UI                        

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
ylim([1 300]);
set(h.trialRaster,'Ydir','reverse');
hold on

% success rate monitor
h.ratePlot = axes(f,'Position',[0.55 0.44 0.3 0.3]);
title('correct rate plot');
xlabel('trial number');     
xlim([1 300]);
ylim([0 1]);
ylabel('correct rate');
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
h.tAirpuff = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.2,...
       'TimerFcn',@(src,event)airpuffEnd);
    %'StartFcn',@(src,event)airpuffStart,...

%% 3 minutes countdown 
%countDown;

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
mLatency = zeros(nTrial,6);   %latency test
%% main loop
start(h.tLickCounter);   %start the lick counter
for trialNum = 1:nTrial
       h.trialNum = trialNum;
    if strcmp(h.tRefractory.Running,'on')
        stop(h.tRefractory);
    end                                      %stop the tRefractory at the beginning of a block
stop(h.tLickCounter)
h.outRWCounterSingle = 0;
h.lickInRWOneTrial = 0;
h.visiOri = h.oriSequence(h.trialNum); %index the visi orientation using trial number
h.targetFlag = (h.visiOri == h.oriTarget); % if target, targetFlag is 1; if not target Flag is 0
if h.targetFlag     %indicator box of target/nontarget trials
   set(h.targetBox,'Visible','on');
elseif ~h.targetFlag
    set(h.nontargetBox,'Visible','on');
end

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
  h.correctRate = (h.hitTrialNumber + h.CRTrialNumber) / trialNum;
  h.data1(trialNum,1) = trialNum;
  h.data1(trialNum,2) = h.resultFlag;   %save the result in the data1 matrix: 1.hit 2.miss 3.FA 4.CR
  h.data1(trialNum,3) = h.lickInRWOneTrial;
  h.data1(trialNum,4) = h.outRWCounterSingle;
  h.data1(trialNum,5) = round(toc(h.trialGlobalTic),2); %time length for individual trials
  h.data1(trialNum,6) = round(h.correctRate,2);
  
  plot(h.ratePlot,trialNum,h.correctRate,'-ok');
  
  if h.targetFlag     %reset target/non-target indicator
   set(h.targetBox,'Visible','off');
 elseif ~h.targetFlag
    set(h.nontargetBox,'Visible','off');
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
disp('---------------------finished!---------------')
totalTime = toc(totalTic);
totalLickTimes = 0;
for i = 1 : numel(h.lickdata)
    totalLickTimes = totalLickTimes + length(h.lickdata{i});
end
PsychPortAudio('Close',h.audioHandle); %close the audio stimulation port
fprintf('>> total time cost:  %s minutes %s seconds \n',num2str(floor((totalTime)/60)),num2str(mod(totalTime,60)));
fprintf('>> total lick times: %s times \n',num2str(h.inRWCounter));
%h.lickdata
%h.data1
save h
sca

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
    while toc(h.postCueTime) <= 1 
    end
    disp('>>post-cue period finished! Visual stimulation starts!')
end

function visiStim(~,~)
  global h
  h.inOrOutRW = 0; 
  updateVbl;
  VSLength = 1; %visual stimulation in the first second
  while h.vbl - h.vblt0 <= VSLength
  Screen('DrawTextures', h.window, h.gabortex, [], [],h.visiOri, [], [], [], [],...
        kPsychDontDoRotation, h.propertiesMat');
  h.vbl = Screen('Flip', h.window, h.vbl + (h.waitframes - 0.5) * h.ifi);
  h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
  end
  updateVbl;
  h.inOrOutRW = 1;
  disp('>>Visual stimulation ended!! Response window starts!!')
end

function responseWindow(~,~)
global h
 h.rwTime = tic;
 h.rwLimit = 4;
 %h.inOrOutRW = 1;  %set as true when RW starts  
  while toc(h.rwTime) <= h.rwLimit    %the RW last for 4 seconds
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
    while toc(h.ITITime) <= h.ITIperiod(h.trialNum) %index the randomized ITI period
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
                   writeDigitalPin(h.a,'D9',1);
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
      stop(h.tLickCounter);
      %disp('timer 1 stopped by timer-2 StartFcn');
      writeDigitalPin(h.a,'D9',0);
end

function refEnd(~,~)
    global h
    start(h.tLickCounter);
    stop(h.tRefractory);
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
 % start(h.tLickCounter);
  % stop(h.tRefractory);
  stop(h.tAirpuff);
end

function seq = randomSequence(n, m)
seq = zeros(n, 1);
seq(1 : m) = 1;
seq(randperm(n)) = seq;
% Up to three consecutive repetitions of the same ID are allowed.
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
        seq(i + 1) = temp; %swap the 4th repetitive element with the end of the sequence

        if seq(i) == seq(i + 1)
            j = j + 1;
        else
            j = 0;
            iSame = 0;
        end
    end
  end
end
