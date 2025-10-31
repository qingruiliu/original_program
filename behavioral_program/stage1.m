%% Stage1 Multi-Session Training Program with UI Monitor
% Modified to include 3 sessions with 100 rewards each
% UI monitoring and camera preview functionality
% Session intervals with black screen countdown

timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor, display the gray color background
global h 
PsychDefaultSetup(2);
Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);
h.screenNumber = max(Screen('Screens'));
h.white = WhiteIndex(h.screenNumber);
h.grey = h.white / 2;
h.black = BlackIndex(h.screenNumber);

[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution); 
h.ifi = Screen('GetFlipInterval',h.window);
h.topPriorityLevel = MaxPriority(h.window);
Priority(h.topPriorityLevel);

%%
%start communication
h.a = arduino("/dev/ttyACM0",'Leonardo','BaudRate',115200);
h.sensorPin = 'D13';
h.waterPumpPin = 'D9';

%%
%display UI, waiting for the initialization
infoUI();

%% Setup UI monitoring interface
setupUI();
Screen('FillRect', h.window, h.black);
Screen('Flip', h.window);

msg = msgbox('Prepare the animal and ready to start experiment.');
waitfor(msg)

%% Initial countdown before starting sessions
initialCountdown();

%% Program variables
h.totalSessions = 3;
h.rewardsPerSession = 100;
h.sessionTimeLimit = 900; % 15 minutes time limit per session (in seconds)
h.sessionIntervalTime = 180; % 3 minutes between sessions
h.allSessionData = {};
h.totalRewards = 0;
h.sessionActive = false; % Flag to track if session is active

%% Define timers
h.tLickCounter = timer('ExecutionMode', 'fixedRate', 'Period', 0.01,...
                             'TimerFcn',@(src,event)pinStatusChanged);
h.tRefractory = timer('BusyMode','error','TasksToExecute',1,'StartDelay',0.1,...
    'StartFcn',@(src,event)refStart,...
    'TimerFcn',@(src,event)refEnd);

%% Main multi-session loop
for sessionNum = 1:h.totalSessions
    fprintf('\n==================== SESSION %d ====================\n', sessionNum);
    h.currentSession = sessionNum;
    h.sessionRewards = 0;
    h.sessionData = zeros(h.rewardsPerSession, 2);
    h.sessionStartTime = tic;
    h.lickFlag = true;
    h.sessionActive = true; % Mark session as active
    
    % Update UI for new session
    updateSessionUI();
    
    % Start lick detection (check if timer is not already running)
    if strcmp(h.tLickCounter.Running, 'off')
        start(h.tLickCounter);
    end
    
    % Session loop - continue until 100 rewards achieved OR time limit reached
    while h.sessionRewards < h.rewardsPerSession && h.sessionActive && toc(h.sessionStartTime) < h.sessionTimeLimit
        pause(0.001); % Small pause to prevent CPU overload

        % Update time display every 10 seconds
        if mod(round(toc(h.sessionStartTime)), 10) == 0
            updateSessionUI();
        end
    end
    
    % Check if session ended due to time limit
    sessionTime = toc(h.sessionStartTime);
    if sessionTime >= h.sessionTimeLimit && h.sessionRewards < h.rewardsPerSession
        fprintf('Session %d ended due to 15-minute time limit! Only %d rewards achieved.\n', ...
            sessionNum, h.sessionRewards);
    end
    
    % Stop lick detection
    if strcmp(h.tLickCounter.Running, 'on')
        stop(h.tLickCounter);
    end
    h.sessionActive = false; % Mark session as inactive
    
    % Store session data
    h.allSessionData{sessionNum} = h.sessionData;
    sessionTime = toc(h.sessionStartTime);
    
    fprintf('Session %d completed! %d rewards in %.1f seconds\n', ...
        sessionNum, h.sessionRewards, sessionTime);
    
    % Update final session stats
    updateSessionUI();
    
    % Inter-session interval (except after last session)
    if sessionNum < h.totalSessions
        fprintf('Starting %d-minute break between sessions...\n', h.sessionIntervalTime/60);
        blackScreenCountdown();
    end
end

%% Clean up and finish
if strcmp(h.tLickCounter.Running, 'on')
    stop(h.tLickCounter);
end
delete(h.tLickCounter);

if strcmp(h.tRefractory.Running, 'on')
    stop(h.tRefractory);
end
delete(h.tRefractory);

fprintf('\n==================== ALL SESSIONS COMPLETED ====================\n');
fprintf('Total rewards across all sessions: %d\n', h.totalRewards);

% Save all data
h.completionTime = datestr(now);
%save h

sca

%% Functions
function infoUI(~,~)
 global h
 prompt = {'mouseID', 'trainStage', 'dayNumber', 'saveDir'};
 dlgtitle = 'mouse information';
 dims = [1 35];
 definput = {'', '', '', '/Users/liuqr/files/MATLAB相关/code ref/test code'};
 h.mouseID = inputdlg(prompt, dlgtitle, dims, definput);
end

function setupUI(~,~)
 global h
 screenSize = get(0,'Screensize'); 
 screenSize(3) = screenSize(3)/2;       
 h.f = figure('Name','Stage1 Multi-Session Monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);
 
 % Session info UI
 uicontrol(h.f,'Style','text','Units','normalized',...
    'Position',[0.05 0.95 0.15 0.04],'String','Current Session','BackgroundColor',[0 1 1],...
    'FontSize',14);
 h.currentSessionUI = uicontrol(h.f,'Style','edit','Units','normalized',...
    'Position',[0.05 0.91 0.15 0.03],'FontSize',14);
 
 uicontrol(h.f,'Style','text','Units','normalized',...
    'Position',[0.22 0.95 0.15 0.04],'String','Session Rewards','BackgroundColor',[0 1 0],...
    'FontSize',14);
 h.sessionRewardsUI = uicontrol(h.f,'Style','edit','Units','normalized',...
    'Position',[0.22 0.91 0.15 0.03],'FontSize',14);
 
 uicontrol(h.f,'Style','text','Units','normalized',...
    'Position',[0.39 0.95 0.15 0.04],'String','Total Rewards','BackgroundColor',[1 1 0],...
    'FontSize',14);
 h.totalRewardsUI = uicontrol(h.f,'Style','edit','Units','normalized',...
    'Position',[0.39 0.91 0.15 0.03],'FontSize',14);
 
 % Session time limit info
 uicontrol(h.f,'Style','text','Units','normalized',...
    'Position',[0.56 0.95 0.18 0.04],'String','Session Time (15min limit)','BackgroundColor',[1 0.7 0.7],...
    'FontSize',14);
 h.sessionTimeUI = uicontrol(h.f,'Style','edit','Units','normalized',...
    'Position',[0.56 0.91 0.18 0.03],'FontSize',14);
 
 % Status indicator
 uicontrol(h.f,'Style','text','Units','normalized',...
    'Position',[0.76 0.95 0.18 0.04],'String','Session Status','BackgroundColor',[0.8 0.8 0.8],...
    'FontSize',14);
 h.sessionStatusUI = uicontrol(h.f,'Style','edit','Units','normalized',...
    'Position',[0.76 0.91 0.18 0.03],'FontSize',14);
 
 % Camera setup (if available)
 try
     % Get camera resolutions
     h.backCam = webcam(1);
     h.frontCam = webcam(2);
     
     backRes = str2double(strsplit(h.backCam.Resolution,'x'));
     frontRes = str2double(strsplit(h.frontCam.Resolution,'x'));
     
     % Calculate aspect ratios
     backAspect = backRes(1) / backRes(2);
     frontAspect = frontRes(1) / frontRes(2);
     
     % Set UI positions based on aspect ratios (adjust for larger camera views)
     backHeight = 0.25;  % Increased height since no plots above
     backWidth = backAspect * backHeight;
     frontHeight = 0.25;  % Increased height since no plots above
     frontWidth = frontAspect * frontHeight;
     
     % Place back camera UI (moved up and made larger)
     h.backCamUI = axes(h.f, 'Position', [0.55 0.3 backWidth backHeight]);
     h.im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', h.backCamUI);
     preview(h.backCam, h.im);
     text(h.backCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 14);
     
     % Place front camera UI (moved up and made larger)
     h.frontCamUI = axes(h.f, 'Position', [0.05 0.3 frontWidth frontHeight]);
     h.im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', h.frontCamUI);
     preview(h.frontCam, h.im2);
     text(h.frontCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 14);
 catch
     fprintf('Warning: Camera setup failed. Continuing without camera preview.\n');
 end
end

function updateSessionUI(~,~)
 global h
 set(h.currentSessionUI, 'String', sprintf('%d/%d', h.currentSession, h.totalSessions));
 set(h.sessionRewardsUI, 'String', sprintf('%d/%d', h.sessionRewards, h.rewardsPerSession));
 set(h.totalRewardsUI, 'String', num2str(h.totalRewards));
 
 % Update session time and status
 if h.sessionActive
     sessionElapsed = toc(h.sessionStartTime);
     timeRemaining = h.sessionTimeLimit - sessionElapsed;
     
     % Format time as MM:SS
     minutes = floor(sessionElapsed / 60);
     seconds = floor(mod(sessionElapsed, 60));
     remainingMinutes = floor(timeRemaining / 60);
     remainingSeconds = floor(mod(timeRemaining, 60));
     
     set(h.sessionTimeUI, 'String', sprintf('%02d:%02d / 15:00', minutes, seconds));
     
     % Update status based on time remaining
     if timeRemaining <= 60
         set(h.sessionStatusUI, 'String', 'TIME WARNING!', 'BackgroundColor', [1 0.3 0.3]);
     elseif timeRemaining <= 300  % 5 minutes
         set(h.sessionStatusUI, 'String', 'TIME ALERT', 'BackgroundColor', [1 0.7 0.3]);
     else
         set(h.sessionStatusUI, 'String', 'ACTIVE', 'BackgroundColor', [0.3 1 0.3]);
     end
 else
     set(h.sessionTimeUI, 'String', '--:-- / 15:00');
     set(h.sessionStatusUI, 'String', 'INACTIVE', 'BackgroundColor', [0.8 0.8 0.8]);
 end
 
 drawnow;
end

function blackScreenCountdown(~,~)
 global h
 % Change screen to black
 Screen('FillRect', h.window, h.black);
 Screen('Flip', h.window);
 
 fprintf('Black screen countdown started...\n');
 
 for remainingSeconds = h.sessionIntervalTime:-1:1
     if mod(remainingSeconds, 30) == 0 || remainingSeconds <= 10
         fprintf('Time remaining: %d seconds\n', remainingSeconds);
     end
     pause(1);
 end
 
 % Return screen to grey
 Screen('FillRect', h.window, h.grey);
 Screen('Flip', h.window);
 
 fprintf('Break finished! Starting next session...\n');
end

function initialCountdown(~,~)
 global h
 % Change screen to black and start initial countdown
 %Screen('FillRect', h.window, h.black);
 %Screen('Flip', h.window);
 
 countdownDuration = 60; % 1 minute countdown
 fprintf('\n==================== INITIAL COUNTDOWN ====================\n');
 fprintf('Please check mouse and lick spout setup!\n');
 fprintf('Sessions will start in 1 minute...\n\n');
 
 for remainingSeconds = countdownDuration:-1:1
     if mod(remainingSeconds, 10) == 0 || remainingSeconds <= 5
         fprintf('Sessions starting in: %d seconds\n', remainingSeconds);
     end
     pause(1);
 end
 
 % Return screen to grey
 Screen('FillRect', h.window, h.grey);
 Screen('Flip', h.window);
 
 fprintf('\nCountdown finished! Starting first session...\n');
end

function pinStatusChanged(~,~)
 global h
 
 % Check if session is active - ignore licks during interval periods
 if ~h.sessionActive
     return;
 end
 
 if strcmp(h.tRefractory.Running,'on')
     return; % Skip if in refractory period
 end
 
 pinValue = readDigitalPin(h.a, h.sensorPin);
 if pinValue == true && h.lickFlag == true
     % Reward delivery
     writeDigitalPin(h.a, h.waterPumpPin, 1);
     h.lickFlag = false;
     h.sessionRewards = h.sessionRewards + 1;
     h.totalRewards = h.totalRewards + 1;
     
     % Record timing data
     sessionTime = toc(h.sessionStartTime);
     h.sessionData(h.sessionRewards, 1) = h.sessionRewards;
     h.sessionData(h.sessionRewards, 2) = sessionTime;
     
     fprintf('Session %d - Reward %d @ %.2f seconds\n', ...
         h.currentSession, h.sessionRewards, sessionTime);
     
     % Update UI
     updateSessionUI();
     
     % Start refractory period
     if strcmp(h.tRefractory.Running, 'off')
         start(h.tRefractory);
     end
 end
end

function refStart(~,~)
 global h
 if strcmp(h.tLickCounter.Running, 'on')
     stop(h.tLickCounter);
 end
 writeDigitalPin(h.a, h.waterPumpPin, 0);
end

function refEnd(~,~)
 global h
 h.lickFlag = true;
 if strcmp(h.tLickCounter.Running, 'off') && h.sessionActive
     start(h.tLickCounter);
 end
 if strcmp(h.tRefractory.Running, 'on')
     stop(h.tRefractory);
 end
end