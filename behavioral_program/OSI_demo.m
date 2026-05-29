% Extended Visual Stimulation Program
% Displays 2 seconds of visual stimulation followed by 4 seconds of grey background
% Orientation changes by 45 degrees each cycle
% Clear the workspace and the screen
sca;
close all;
clear;
 
% Setup PTB with some default values
PsychDefaultSetup(2);

% initial setting of PTB-3
Screen('Preference', 'ConserveVRAM', 4096); 
Screen('Preference', 'VBLTimestampingMode', 4);  
Screen('Preference', 'SkipSyncTests', 0); 
Screen('Preference', 'VisualDebugLevel', 0); 
Screen('Preference', 'SuppressAllWarnings', 1); 

Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);

% Set the screen number to the external secondary monitor if there is one connected
screenNumber = max(Screen('Screens'));
 
% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;

% Open the screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey,...
    [], 32, 2, [], [], kPsychNeedRetinaResolution);

% Get the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);

% Set maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
for i = 1:5
    Screen('Flip', window);
end

%--------------------
% Gabor information
%--------------------

% Calculate diagonal length to ensure full screen coverage when rotated
screenWidth = windowRect(3);
screenHeight = windowRect(4);
diagonal = sqrt(screenWidth^2 + screenHeight^2);

% Use diagonal length for both width and height to ensure full coverage
width = diagonal;
height = diagonal;
  
% Sigma of Gaussian - adjust based on screen width for consistency
sigma = screenWidth / 7;

% Parameters
contrast = 1;
aspectRatio = 1.0;
phase = 0;

% Spatial Frequency (Cycles Per Pixel)
numCycles = 4;
freq = numCycles / screenWidth;

% Build a procedural gabor texture
backgroundOffset = [0.5 0.5 0.5 0.0];
disableNorm = 1;
preContrastMultiplier = 0.5;
gabortex = CreateProceduralSineGrating(window, width, height, backgroundOffset,...
    [], preContrastMultiplier);

%--------------------
% Experimental parameters
%--------------------

% Timing parameters
stimDuration = 2;    % 2 seconds of visual stimulation
greyDuration = 4;    % 4 seconds of grey background
waitframes = 1;
phasePerFrame = 4 * pi;

% Orientation parameters - 45 degree increments (8 total orientations)
orientations = 0:45:315;  % [0, 45, 90, 135, 180, 225, 270, 315] degrees
numOrientations = length(orientations);
currentOrientationIndex = 1;

% Trial counter
trialNumber = 1;

%% monitor UI
    screenSize = get(0,'Screensize');
    screenSize(3) = screenSize(3)/2; % 使用半屏宽度
    monitorFig = figure('Name', 'Passive Viewing 实验监控', ...
        'Position', screenSize, 'Color', [0.95 0.95 0.95]);
    
    % 摄像头面板
    cameraPanel = uipanel('Parent', monitorFig, 'Title', '摄像头监控', ...
        'Position', [0.02 0.45 0.96 0.38], 'FontSize', 12, 'FontWeight', 'bold');
    
    try
        % 初始化摄像头
        cam1 = webcam(1);
        cam2 = webcam(2);
        cameraAvailable = true;
        
        % 摄像头1 (左侧)
        cam1Axes = axes('Parent', cameraPanel, 'Position', [0.05 0.1 0.4 0.8]);
        cam1Image = image(cam1Axes, snapshot(cam1));
        title(cam1Axes, '摄像头1 - 后方监控');
        axis(cam1Axes, 'off');
        
        % 摄像头2 (右侧)  
        cam2Axes = axes('Parent', cameraPanel, 'Position', [0.55 0.1 0.4 0.8]);
        cam2Image = image(cam2Axes, snapshot(cam2));
        title(cam2Axes, '摄像头2 - 前方监控');
        axis(cam2Axes, 'off');
        
    catch
        h.cameraAvailable = false;
        h.noCamAxes = axes('Parent', h.cameraPanel, 'Position', [0.1 0.1 0.8 0.8]);
        text(h.noCamAxes, 0.5, 0.5, '摄像头未连接或不可用', ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'r');
        axis(h.noCamAxes, [0 1 0 1]);
        axis(h.noCamAxes, 'off');
        warning('摄像头初始化失败，继续进行实验但无摄像头监控。');
    end

%% Display initial information
fprintf('Extended Visual Stimulation Program\n');
fprintf('Orientations: %s degrees\n', mat2str(orientations));
fprintf('Stimulation duration: %d seconds\n', stimDuration);
fprintf('Grey background duration: %d seconds\n', greyDuration);
msgb = msgbox('Press any key to stop the program\n');
waitfor(msgb);
fprintf('==========================================\n');

%------------------------------------------
% Main experimental loop
%------------------------------------------

% Initial flip to start timing
vbl = Screen('Flip', window);
startTime = vbl;

while ~KbCheck
    
    % Get current orientation
    currentOrientation = orientations(currentOrientationIndex);
    
    % Display trial information
    fprintf('Trial %d: Orientation = %d degrees\n', trialNumber, currentOrientation);
    
    %% Visual Stimulation Phase (2 seconds)
    fprintf('  Visual stimulation phase started...\n');
    
    % Reset phase for each trial
    currentPhase = phase;
    
    % Create properties matrix for current orientation
    propertiesMat = [currentPhase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
    
    % Record start time for this stimulation phase
    stimStartTime = vbl;
    
    % Visual stimulation loop
    while (vbl - stimStartTime) < stimDuration
        % Draw the Gabor with current orientation
        Screen('DrawTextures', window, gabortex, [], [], currentOrientation, [], [], [], [],...
            [], propertiesMat');
        
        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % Update the phase for animation
        propertiesMat(1) = propertiesMat(1) + phasePerFrame;
        
        % Check for key press to exit
        if KbCheck
            break;
        end
    end
    
    %% Grey Background Phase (4 seconds)
    fprintf('  Grey background phase started...\n');
    
    % Clear screen to grey background
    Screen('FillRect', window, grey);
    vbl = Screen('Flip', window);
    greyStartTime = vbl;
    
    % Grey background loop
    while (vbl - greyStartTime) < greyDuration
        % Keep grey background
        Screen('FillRect', window, grey);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % Check for key press to exit
        if KbCheck
            break;
        end
    end
    
    % Check for key press to exit (final check)
    if KbCheck
        break;
    end
    
    % Move to next orientation
    currentOrientationIndex = currentOrientationIndex + 1;
    if currentOrientationIndex > numOrientations
        currentOrientationIndex = 1;  % Reset to first orientation
    end
    
    % Increment trial number
    trialNumber = trialNumber + 1;
    
    fprintf('  Trial %d completed. Total time: %.2f seconds\n', trialNumber-1, vbl - stimStartTime);
    fprintf('------------------------------------------\n');
    
end

% Calculate total experimental time
totalTime = vbl - startTime;
fprintf('\nExperiment completed!\n');
fprintf('Total trials: %d\n', trialNumber-1);
fprintf('Total time: %.2f seconds (%.2f minutes)\n', totalTime, totalTime/60);

% Clear screen
sca;