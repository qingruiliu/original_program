% New Visual Stimulation Program with Noise and AB Pattern
% Trial sequence: Random noise (2s) -> Grey (2s) -> ABAB pattern (A=120°, B=90°)
% A and B are 2 seconds each followed by 4 seconds grey interval
% 100 trials total
% Clear the workspace and the screen
sca;
close all;
clear;
 
% Setup PTB with some default values
PsychDefaultSetup(2); 

% Initial setting of PTB-3
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
black = BlackIndex(screenNumber);
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
% Screen parameters
%--------------------
screenWidth = windowRect(3);
screenHeight = windowRect(4);
diagonal = sqrt(screenWidth^2 + screenHeight^2);

%--------------------
% Noise stimulus parameters
%--------------------
noiseFrameDuration = 0.5;  % Each noise frame lasts 0.5 seconds
noiseTotalDuration = 2;     % Total noise duration: 2 seconds
numNoiseFrames = noiseTotalDuration / noiseFrameDuration;  % 4 frames
blockSize = 80;  % Size of each noise block in pixels
numBlackBlocks = 10;  % Number of black blocks per frame
numWhiteBlocks = 10;  % Number of white blocks per frame

%--------------------
% Gabor information for oriented stimuli
%--------------------

% Use diagonal length for both width and height to ensure full coverage
width = diagonal;
height = diagonal;
  
% Sigma of Gaussian
sigma = screenWidth / 7;

% Parameters
contrast = 1;
aspectRatio = 1.0;
phase = 0;

% Spatial Frequency (Cycles Per Pixel)
numCycles = 4;
freq = numCycles / screenWidth;

% Build a procedural gabor texture for oriented stimuli
backgroundOffset = [0.5 0.5 0.5 0.0];
disableNorm = 1;
preContrastMultiplier = 0.5;
gabortex = CreateProceduralSineGrating(window, width, height, backgroundOffset,...
    [], preContrastMultiplier);

%--------------------
% Experimental parameters
%--------------------

% Timing parameters
noiseDuration = 2;           % 2 seconds of random noise
greyInitialDuration = 2;     % 2 seconds of grey after noise
stimDuration = 2;            % 2 seconds for each oriented stimulus (A or B)
greyIntervalDuration = 4;    % 4 seconds grey interval after each stimulus
waitframes = 1;
phasePerFrame = 4 * pi;

% Orientation parameters
orientationA = 120;  % Stimulus A: 120 degrees
orientationB = 90;   % Stimulus B: 90 degrees

% Trial parameters
maxTrials = 100;  % 100 trials total

% Display initial information
fprintf('New Visual Stimulation Program - Noise + AB Pattern\n');
fprintf('==========================================\n');
fprintf('Trial structure:\n');
fprintf('  1. Random noise blocks (2s, 4 frames @ 0.5s each, 10 black + 10 white blocks)\n');
fprintf('  2. Grey background (2s)\n');
fprintf('  3. Stimulus A (120°, 2s) + Grey interval (4s)\n');
fprintf('  4. Stimulus B (90°, 2s) + Grey interval (4s)\n');
fprintf('  5. Stimulus A (120°, 2s) + Grey interval (4s)\n');
fprintf('  6. Stimulus B (90°, 2s) + Grey interval (4s)\n');
fprintf('Total trials: %d\n', maxTrials);
fprintf('Total experiment duration: ~%.0f seconds\n', maxTrials * (2 + 2 + 4*(2+4)));
fprintf('==========================================\n');

%% Create GUI for monitoring
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','Trial Monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);   

% Add trial information display
infoPanel = uipanel(f, 'Position', [0.05 0.45 0.9 0.5], 'Title', 'Trial Information', 'FontSize', 12);
trialText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Waiting to start...', ...
    'Units', 'normalized', 'Position', [0.05 0.7 0.9 0.25], ...
    'FontSize', 14, 'HorizontalAlignment', 'left', 'BackgroundColor', [0.95 0.95 0.95]);
statusText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Status: Ready', ...
    'Units', 'normalized', 'Position', [0.05 0.4 0.9 0.25], ...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'BackgroundColor', [0.95 0.95 0.95]);
timeText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Elapsed time: 0.0 s', ...
    'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.25], ...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'BackgroundColor', [0.95 0.95 0.95]);

% Initialize camera variables
cameraInitialized = false;

% Open cameras with error handling
try
    backCam = webcam(1);
    frontCam = webcam(2);

    backRes = str2double(strsplit(backCam.Resolution,'x'));
    frontRes = str2double(strsplit(frontCam.Resolution,'x'));

    % Calculate aspect ratios
    backAspect = backRes(1) / backRes(2);
    frontAspect = frontRes(1) / frontRes(2);

    % Set UI positions based on aspect ratios
    backHeight = 0.28;
    backWidth = backAspect * backHeight;
    frontHeight = 0.28;
    frontWidth = frontAspect * frontHeight;

    % Place back camera UI
    backCamUI = axes(f, 'Position', [0.45 0.1 backWidth backHeight]);
    im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', backCamUI);
    preview(backCam, im);
    text(backCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);

    % Place front camera UI
    frontCamUI = axes(f, 'Position', [0.05 0.1 frontWidth frontHeight]);
    im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', frontCamUI);
    preview(frontCam, im2);
    text(frontCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    
    cameraInitialized = true;
    fprintf('Camera initialization successful.\n');
    
    % Lower priority slightly to allow camera operation
    Priority(0);
catch ME
    warning('Camera initialization failed: %s. Continuing without cameras.', ME.message);
    cameraInitialized = false;
end

%% Confirmation and countdown before experiment
msg = msgbox('Ready to start! Press OK to begin countdown.');
waitfor(msg)
fprintf('Starting 30-second countdown...\n');

% 30-second countdown
for countdown = 30:-1:1
    fprintf('Time remaining: %d seconds\n', countdown);
    pause(1);
end

% Clear screen
Screen('FillRect', window, grey);
Screen('Flip', window);

fprintf('Countdown completed. Experiment starting now!\n');
fprintf('Press any key to stop the program\n');
fprintf('==========================================\n');

%% Initialize timing recording
tic;  % Start master timer
masterStartTime = toc;

% Initialize arrays to store timestamps and event information
trialStartTimes = [];
noiseStartTimes = [];
greyInitialStartTimes = [];
stimAStartTimes = [];
stimBStartTimes = [];
greyIntervalStartTimes = [];
eventLog = {};  % Cell array to log all events

%% Main experimental loop
try
    % Set high priority for visual stimulation timing precision
    Priority(topPriorityLevel);
    
    % Initial flip to start timing
    vbl = Screen('Flip', window);
    experimentStartTime = vbl;

    for trialNum = 1:maxTrials
        
        % Update GUI
        set(trialText, 'String', sprintf('Trial: %d / %d', trialNum, maxTrials));
        drawnow limitrate;
        
        % Record trial start
        trialStartTime = toc;
        trialStartTimes = [trialStartTimes; trialStartTime];
        fprintf('\n========== TRIAL %d / %d ==========\n', trialNum, maxTrials);
        fprintf('[%.3f s] Trial %d started\n', trialStartTime, trialNum);
        eventLog{end+1} = sprintf('%.3f,Trial %d Start', trialStartTime, trialNum);
        
        %% Phase 1: Random Noise (2 seconds, 4 frames of 0.5s each)
        set(statusText, 'String', 'Status: Random Noise');
        drawnow limitrate;
        
        noiseStartTime = toc;
        noiseStartTimes = [noiseStartTimes; noiseStartTime];
        fprintf('[%.3f s] Random noise phase started\n', noiseStartTime);
        eventLog{end+1} = sprintf('%.3f,Noise Start', noiseStartTime);
        
        for noiseFrame = 1:numNoiseFrames
            % Start with grey background
            Screen('FillRect', window, grey);
            
            % Store positions of black and white blocks separately
            blackBlockRects = [];
            whiteBlockRects = [];
            
            % Generate random positions for black blocks (10 blocks, non-overlapping with each other)
            for i = 1:numBlackBlocks
                overlap = true;
                attempts = 0;
                maxAttempts = 100;  % Prevent infinite loop
                
                while overlap && attempts < maxAttempts
                    % Random position ensuring block stays within screen
                    x = randi([0, screenWidth - blockSize]);
                    y = randi([0, screenHeight - blockSize]);
                    blockRect = [x, y, x + blockSize, y + blockSize];
                    
                    % Check if this block overlaps with existing black blocks
                    overlap = false;
                    for j = 1:size(blackBlockRects, 1)
                        existingRect = blackBlockRects(j, :);
                        % Check rectangle overlap
                        if ~(blockRect(3) <= existingRect(1) || blockRect(1) >= existingRect(3) || ...
                             blockRect(4) <= existingRect(2) || blockRect(2) >= existingRect(4))
                            overlap = true;
                            break;
                        end
                    end
                    
                    attempts = attempts + 1;
                end
                
                % Store and draw black block
                blackBlockRects = [blackBlockRects; blockRect];
                Screen('FillRect', window, black, blockRect);
            end
            
            % Generate random positions for white blocks (10 blocks, non-overlapping with each other)
            for i = 1:numWhiteBlocks
                overlap = true;
                attempts = 0;
                maxAttempts = 100;  % Prevent infinite loop
                
                while overlap && attempts < maxAttempts
                    % Random position ensuring block stays within screen
                    x = randi([0, screenWidth - blockSize]);
                    y = randi([0, screenHeight - blockSize]);
                    blockRect = [x, y, x + blockSize, y + blockSize];
                    
                    % Check if this block overlaps with existing white blocks
                    overlap = false;
                    for j = 1:size(whiteBlockRects, 1)
                        existingRect = whiteBlockRects(j, :);
                        % Check rectangle overlap
                        if ~(blockRect(3) <= existingRect(1) || blockRect(1) >= existingRect(3) || ...
                             blockRect(4) <= existingRect(2) || blockRect(2) >= existingRect(4))
                            overlap = true;
                            break;
                        end
                    end
                    
                    attempts = attempts + 1;
                end
                
                % Store and draw white block
                whiteBlockRects = [whiteBlockRects; blockRect];
                Screen('FillRect', window, white, blockRect);
            end
            
            % Flip to show the noise pattern and hold it for 0.5 seconds
            vbl = Screen('Flip', window);
            frameStartTime = vbl;
            
            % Hold this frame for 0.5 seconds
            while (vbl - frameStartTime) < noiseFrameDuration
                % Check for key press to exit
                if KbCheck
                    break;
                end
                
                % Small pause to reduce CPU usage while waiting
                WaitSecs(0.01);
                vbl = GetSecs();
                
                % Allow GUI updates periodically
                if mod(round((vbl - frameStartTime) * 100), 10) == 0
                    drawnow limitrate;
                end
            end
            
            % Check for key press to exit
            if KbCheck
                break;
            end
        end
        
        % Check for key press to exit
        if KbCheck
            break;
        end
        
        %% Phase 2: Grey Background (2 seconds)
        set(statusText, 'String', 'Status: Initial Grey Background');
        drawnow limitrate;
        
        greyInitialStartTime = toc;
        greyInitialStartTimes = [greyInitialStartTimes; greyInitialStartTime];
        fprintf('[%.3f s] Initial grey background phase started\n', greyInitialStartTime);
        eventLog{end+1} = sprintf('%.3f,Initial Grey Start', greyInitialStartTime);
        
        Screen('FillRect', window, grey);
        vbl = Screen('Flip', window);
        greyStartTime = vbl;
        
        frameCount = 0;
        while (vbl - greyStartTime) < greyInitialDuration
            frameCount = frameCount + 1;
            
            Screen('FillRect', window, grey);
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
            % Allow GUI updates
            if mod(frameCount, 5) == 0
                set(timeText, 'String', sprintf('Elapsed time: %.0f s', toc));
                drawnow limitrate;
            end
            
            % Check for key press to exit
            if KbCheck
                break;
            end
        end
        
        % Check for key press to exit
        if KbCheck
            break;
        end
        
        %% Phase 3-6: ABAB pattern (A=120°, B=90°)
        % Each stimulus is 2s followed by 4s grey interval
        % Pattern: A -> grey -> B -> grey -> A -> grey -> B -> grey
        
        stimSequence = [orientationA, orientationB, orientationA, orientationB];
        stimNames = {'A (120°)', 'B (90°)', 'A (120°)', 'B (90°)'};
        
        for stimIdx = 1:4
            currentOrientation = stimSequence(stimIdx);
            stimName = stimNames{stimIdx};
            
            %% Oriented Stimulus (2 seconds)
            set(statusText, 'String', sprintf('Status: Stimulus %s', stimName));
            drawnow limitrate;
            
            stimStartTime = toc;
            if currentOrientation == orientationA
                stimAStartTimes = [stimAStartTimes; stimStartTime];
            else
                stimBStartTimes = [stimBStartTimes; stimStartTime];
            end
            
            fprintf('[%.3f s] Stimulus %s started\n', stimStartTime, stimName);
            eventLog{end+1} = sprintf('%.3f,Stimulus %s Start', stimStartTime, stimName);
            
            % Reset phase
            currentPhase = phase;
            
            % Create properties matrix for current orientation
            propertiesMat = [currentPhase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
            
            % Stimulus display loop
            stimPhaseStartTime = vbl;
            frameCount = 0;
            while (vbl - stimPhaseStartTime) < stimDuration
                frameCount = frameCount + 1;
                
                % Draw the Gabor with current orientation
                Screen('DrawTextures', window, gabortex, [], [], currentOrientation, [], [], [], [],...
                    [], propertiesMat');
                
                % Flip to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                % Update the phase for animation
                propertiesMat(1) = propertiesMat(1) + phasePerFrame;
                
                % Allow GUI updates
                if mod(frameCount, 10) == 0
                    drawnow limitrate;
                end
                
                % Check for key press to exit
                if KbCheck
                    break;
                end
            end
            
            % Check for key press to exit
            if KbCheck
                break;
            end
            
            %% Grey Interval (4 seconds)
            set(statusText, 'String', 'Status: Grey Interval');
            drawnow limitrate;
            
            greyIntervalStartTime = toc;
            greyIntervalStartTimes = [greyIntervalStartTimes; greyIntervalStartTime];
            fprintf('[%.3f s] Grey interval after stimulus %s started\n', greyIntervalStartTime, stimName);
            eventLog{end+1} = sprintf('%.3f,Grey Interval Start (after %s)', greyIntervalStartTime, stimName);
            
            Screen('FillRect', window, grey);
            vbl = Screen('Flip', window);
            greyIntervalPhaseStartTime = vbl;
            
            frameCount = 0;
            while (vbl - greyIntervalPhaseStartTime) < greyIntervalDuration
                frameCount = frameCount + 1;
                
                Screen('FillRect', window, grey);
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                
                % Allow GUI updates
                if mod(frameCount, 5) == 0
                    set(timeText, 'String', sprintf('Elapsed time: %.1f s', toc));
                    drawnow limitrate;
                end
                
                % Check for key press to exit
                if KbCheck
                    break;
                end
            end
            
            % Check for key press to exit
            if KbCheck
                break;
            end
        end
        
        % Check for key press to exit
        if KbCheck
            break;
        end
        
        % Trial completed
        trialEndTime = toc;
        fprintf('[%.3f s] Trial %d completed (Duration: %.2f s)\n', ...
            trialEndTime, trialNum, trialEndTime - trialStartTime);
        eventLog{end+1} = sprintf('%.3f,Trial %d End', trialEndTime, trialNum);
        
    end  % End of trial loop
    
catch ME
    % Error handling
    fprintf('\nERROR occurred during experiment:\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('In file: %s\n', ME.stack(1).file);
    fprintf('At line: %d\n', ME.stack(1).line);
    
    % Display error information
    Priority(0);
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end

%% Experiment completed - cleanup and save data
Priority(0);

fprintf('\n==========================================\n');
fprintf('Experiment completed successfully!\n');
fprintf('Total trials completed: %d\n', trialNum);
fprintf('Total experiment duration: %.0f seconds\n', toc);
fprintf('==========================================\n');

% Save timing data
fprintf('Saving experiment data...\n');
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
filename = sprintf('noise_AB_experiment_%s.mat', currentDateTime);

% Save all timing and event data
save(filename, 'trialStartTimes', 'noiseStartTimes', 'greyInitialStartTimes', ...
    'stimAStartTimes', 'stimBStartTimes', 'greyIntervalStartTimes', ...
    'eventLog', 'maxTrials', 'orientationA', 'orientationB', ...
    'noiseDuration', 'greyInitialDuration', 'stimDuration', 'greyIntervalDuration', ...
    'masterStartTime');

fprintf('Data saved to: %s\n', filename);

% Save event log to CSV file
csvFilename = sprintf('noise_AB_experiment_log_%s.csv', currentDateTime);
fid = fopen(csvFilename, 'w');
fprintf(fid, 'Time(s),Event\n');
for i = 1:length(eventLog)
    fprintf(fid, '%s\n', eventLog{i});
end
fclose(fid);
fprintf('Event log saved to: %s\n', csvFilename);

% Close cameras if initialized
if cameraInitialized
    try
        clear backCam frontCam;
        fprintf('Cameras closed successfully.\n');
    catch
        warning('Failed to close cameras properly.');
    end
end

% Close the figure
try
    close(f);
catch
    % Figure might already be closed
end

% Close the screen
sca;

fprintf('\nProgram completed. All resources cleaned up.\n');
fprintf('Thank you for using the experiment program!\n');
