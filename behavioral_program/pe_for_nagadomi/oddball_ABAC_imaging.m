% Prediction Error Detection Program with ABAB/ABAC Pattern
% Trial sequence: Random noise (2s) -> Grey (2s) -> ABAB or ABAC pattern
% A=120°, B=90°, C=240° (oddball stimulus)
% First 60 trials: 20% probability of ABAC (no C in first 5 trials, max 2 consecutive)
% Last 40 trials: All ABAC pattern
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
stimDuration = 2;            % 2 seconds for each oriented stimulus (A, B, or C)
greyIntervalDuration = 4;    % 4 seconds grey interval after each stimulus
waitframes = 1;
phasePerFrame = 4 * pi;

% Orientation parameters
orientationA = 120;  % Stimulus A: 120 degrees
orientationB = 90;   % Stimulus B: 90 degrees
orientationC = 240;  % Stimulus C: 240 degrees (oddball)

% Trial parameters
maxTrials = 100;  % 100 trials total
habituationTrials = 60;  % First 60 trials with 20% oddball
testTrials = 40;  % Last 40 trials all with oddball

%--------------------
% Generate trial sequence
%--------------------

% Initialize trial types array (0 = ABAB, 1 = ABAC)
trialTypes = zeros(maxTrials, 1);

% First 3 trials: all ABAB (no C)
% (already initialized to 0)

% Trials 4-60: 20% of first 60 trials = 12 ABAC trials with constraint of max 2 consecutive
targetOddballs = round(60 * 0.2);  % 12 oddball trials in trials 1-60
oddballsPlaced = 0;
consecutiveOddballs = 0;

% Shuffle trials 4-60 to place oddballs
availableTrials = 4:60;
availableTrials = availableTrials(randperm(length(availableTrials)));

for i = 1:length(availableTrials)
    trialIdx = availableTrials(i);
    
    if oddballsPlaced < targetOddballs
        % Check if we can place an oddball here
        if consecutiveOddballs < 2
            trialTypes(trialIdx) = 1;  % ABAC
            oddballsPlaced = oddballsPlaced + 1;
            consecutiveOddballs = consecutiveOddballs + 1;
        else
            % Must place ABAB to break consecutive streak
            consecutiveOddballs = 0;
        end
    else
        % Reset consecutive counter for non-oddball
        if trialTypes(trialIdx) == 0
            consecutiveOddballs = 0;
        end
    end
end

% Trials 61-100: all ABAC
trialTypes(61:100) = 1;

% Display initial information
fprintf('Prediction Error Detection Program - ABAB/ABAC Pattern\n');
fprintf('==========================================\n');
fprintf('Trial structure:\n');
fprintf('  1. Random noise blocks (2s, 4 frames @ 0.5s each, 10 black + 10 white blocks)\n');
fprintf('  2. Grey background (2s)\n');
fprintf('  3. Stimulus A (120°, 2s) + Grey interval (4s)\n');
fprintf('  4. Stimulus B (90°, 2s) + Grey interval (4s)\n');
fprintf('  5. Stimulus A (120°, 2s) + Grey interval (4s)\n');
fprintf('  6. Stimulus B or C (90° or 240°, 2s) + Grey interval (4s)\n');
fprintf('==========================================\n');
fprintf('Trial sequence:\n');
fprintf('  Trials 1-3: All ABAB (no oddball)\n');
fprintf('  Trials 4-60: 20%% ABAC (oddball C at 240°), max 2 consecutive\n');
fprintf('  Trials 61-100: All ABAC (oddball C at 240°)\n');
fprintf('Total trials: %d\n', maxTrials);
fprintf('Oddball trials in 1-60: %d / 60 (%.1f%%)\n', sum(trialTypes(1:60)), sum(trialTypes(1:60))/60*100);
fprintf('Total experiment duration: ~%.0f seconds\n', maxTrials * (2 + 2 + 4*(2+4)));
fprintf('==========================================\n');

%% Create GUI for monitoring
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','Prediction Error Monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);   

% Add trial information display
infoPanel = uipanel(f, 'Position', [0.05 0.45 0.9 0.5], 'Title', 'Trial Information', 'FontSize', 12);
trialText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Waiting to start...', ...
    'Units', 'normalized', 'Position', [0.05 0.7 0.9 0.25], ...
    'FontSize', 14, 'HorizontalAlignment', 'left', 'BackgroundColor', [0.95 0.95 0.95]);
statusText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Status: Ready', ...
    'Units', 'normalized', 'Position', [0.05 0.4 0.9 0.25], ...
    'FontSize', 12, 'HorizontalAlignment', 'left', 'BackgroundColor', [0.95 0.95 0.95]);
timeText = uicontrol(infoPanel, 'Style', 'text', 'String', 'Elapsed time: 0 s', ...
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
stimCStartTimes = [];
greyIntervalStartTimes = [];
eventLog = {};  % Cell array to log all events
actualTrialTypes = [];  % Record actual trial types

%% Main experimental loop
try
    % Set high priority for visual stimulation timing precision
    Priority(topPriorityLevel);
    
    % Initial flip to start timing
    vbl = Screen('Flip', window);
    experimentStartTime = vbl;

    for trialNum = 1:maxTrials
        
        % Determine trial type for this trial
        isOddballTrial = trialTypes(trialNum);
        trialTypeStr = 'ABAB';
        if isOddballTrial
            trialTypeStr = 'ABAC (ODDBALL)';
        end
        
        % Update GUI
        set(trialText, 'String', sprintf('Trial: %d / %d - Type: %s', trialNum, maxTrials, trialTypeStr));
        drawnow limitrate;
        
        % Record trial start
        trialStartTime = toc;
        trialStartTimes = [trialStartTimes; trialStartTime];
        actualTrialTypes = [actualTrialTypes; isOddballTrial];
        fprintf('\n========== TRIAL %d / %d - %s ==========\n', trialNum, maxTrials, trialTypeStr);
        fprintf('[%.3f s] Trial %d started\n', trialStartTime, trialNum);
        eventLog{end+1} = sprintf('%.3f,Trial %d Start,%s', trialStartTime, trialNum, trialTypeStr);
        
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
        
        %% Phase 3-6: ABAB or ABAC pattern
        % Each stimulus is 2s followed by 4s grey interval
        % Pattern: A -> grey -> B -> grey -> A -> grey -> B/C -> grey
        
        % Determine the sequence based on trial type
        if isOddballTrial
            stimSequence = [orientationA, orientationB, orientationA, orientationC];
            stimNames = {'A (120°)', 'B (90°)', 'A (120°)', 'C (240°) ODDBALL'};
        else
            stimSequence = [orientationA, orientationB, orientationA, orientationB];
            stimNames = {'A (120°)', 'B (90°)', 'A (120°)', 'B (90°)'};
        end
        
        for stimIdx = 1:4
            currentOrientation = stimSequence(stimIdx);
            stimName = stimNames{stimIdx};
            
            %% Oriented Stimulus (2 seconds)
            set(statusText, 'String', sprintf('Status: Stimulus %s', stimName));
            drawnow limitrate;
            
            stimStartTime = toc;
            if currentOrientation == orientationA
                stimAStartTimes = [stimAStartTimes; stimStartTime];
            elseif currentOrientation == orientationB
                stimBStartTimes = [stimBStartTimes; stimStartTime];
            else  % orientationC
                stimCStartTimes = [stimCStartTimes; stimStartTime];
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
fprintf('Total oddball (ABAC) trials: %d\n', sum(actualTrialTypes));
fprintf('Total standard (ABAB) trials: %d\n', sum(~actualTrialTypes));
fprintf('Total experiment duration: %.0f seconds\n', toc);
fprintf('==========================================\n');

% Save timing data
fprintf('Saving experiment data...\n');
currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
filename = sprintf('prediction_error_experiment_%s.mat', currentDateTime);

% Save all timing and event data
save(filename, 'trialStartTimes', 'noiseStartTimes', 'greyInitialStartTimes', ...
    'stimAStartTimes', 'stimBStartTimes', 'stimCStartTimes', 'greyIntervalStartTimes', ...
    'eventLog', 'maxTrials', 'orientationA', 'orientationB', 'orientationC', ...
    'noiseDuration', 'greyInitialDuration', 'stimDuration', 'greyIntervalDuration', ...
    'trialTypes', 'actualTrialTypes', 'masterStartTime');

fprintf('Data saved to: %s\n', filename);

% Save event log to CSV file
csvFilename = sprintf('prediction_error_experiment_log_%s.csv', currentDateTime);
fid = fopen(csvFilename, 'w');
fprintf(fid, 'Time(s),Event\n');
for i = 1:length(eventLog)
    fprintf(fid, '%s\n', eventLog{i});
end
fclose(fid);
fprintf('Event log saved to: %s\n', csvFilename);

% Save trial sequence to CSV file
trialSeqFilename = sprintf('prediction_error_trial_sequence_%s.csv', currentDateTime);
fid = fopen(trialSeqFilename, 'w');
fprintf(fid, 'TrialNumber,TrialType,IsOddball\n');
for i = 1:length(actualTrialTypes)
    if actualTrialTypes(i)
        fprintf(fid, '%d,ABAC,1\n', i);
    else
        fprintf(fid, '%d,ABAB,0\n', i);
    end
end
fclose(fid);
fprintf('Trial sequence saved to: %s\n', trialSeqFilename);

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
