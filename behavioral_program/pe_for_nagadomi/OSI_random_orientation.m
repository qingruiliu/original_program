% Extended Visual Stimulation Program with Randomized Orientation Order
% Displays 2 seconds of visual stimulation followed by 4 seconds of grey background
% Orientation order is randomized within each cycle
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
greyDuration = 4;    % 4 seconds of grey background between stimuli
cycleBreakDuration = 6;  % 6 seconds break between complete orientation cycles
waitframes = 1;
phasePerFrame = 4 * pi;

% Orientation parameters - 30 degree increments (12 total orientations)
baseOrientations = 0:30:330;  % [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330] degrees
numOrientations = length(baseOrientations);

% Cycle and trial parameters
maxCycles = 20;  % 总共20个循环
currentCycle = 1;  % 当前循环数
maxTrials = maxCycles * numOrientations;  % 总试次数 = 20 × 12 = 240

% Trial counter
trialNumber = 1;

% Initialize orientation sequence for first cycle (randomized)
currentCycleOrientations = baseOrientations(randperm(numOrientations));
currentOrientationIndex = 1;

% Display initial information
fprintf('Extended Visual Stimulation Program - Randomized Orientation Order\n');
fprintf('Base Orientations: %s degrees (30-degree increments)\n', mat2str(baseOrientations));
fprintf('Orientation order is RANDOMIZED within each cycle\n');
fprintf('Total cycles: %d (12 stimuli per cycle)\n', maxCycles);
fprintf('Total trials: %d\n', maxTrials);
fprintf('Stimulation duration: %d seconds\n', stimDuration);
fprintf('Grey background duration: %d seconds\n', greyDuration);
fprintf('Cycle break duration: %d seconds\n', cycleBreakDuration);
fprintf('==========================================\n');
fprintf('Cycle 1 orientation order: %s\n', mat2str(currentCycleOrientations));
fprintf('==========================================\n');

%% create GUI
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);   

% Initialize camera variables
cameraInitialized = false;

% Open cameras with error handling - 使用stage3的双摄像头功能
try
    backCam = webcam(1);
    frontCam = webcam(2);

    backRes = str2double(strsplit(backCam.Resolution,'x'));
    frontRes = str2double(strsplit(frontCam.Resolution,'x'));

    % Calculate aspect ratios
    backAspect = backRes(1) / backRes(2);
    frontAspect = frontRes(1) / frontRes(2);

    % Set UI positions based on aspect ratios (normalized units) - 调整位置以适应新布局
    backHeight = 0.28;
    backWidth = backAspect * backHeight;
    frontHeight = 0.28;
    frontWidth = frontAspect * frontHeight;

    % Place back camera UI - 调整位置
    backCamUI = axes(f, 'Position', [0.45 0.1 backWidth backHeight]);
    im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', backCamUI);
    preview(backCam, im);
    text(backCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);  % Add label below the image

    % Place front camera UI - 调整位置
    frontCamUI = axes(f, 'Position', [0.05 0.1 frontWidth frontHeight]);
    title(frontCamUI, 'Back Camera');
    im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', frontCamUI);
    preview(frontCam, im2);
    text(frontCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);  % Add label below the image
    
    cameraInitialized = true;
    fprintf('Camera initialization successful.\n');
    
    % Lower priority slightly to allow camera operation
    Priority(0);  % Reset to normal priority after camera setup
catch ME
    warning('Camera initialization failed: %s. Continuing without cameras.', ME.message);
    cameraInitialized = false;
end
%%
%------------------------------------------
% Confirmation and countdown before experiment
%------------------------------------------
msg =msgbox('start!');
waitfor(msg)
fprintf('Starting 30-second countdown...\n');

% 30-second countdown - 只在命令行显示
for countdown = 30:-1:1
 
     fprintf('Time remaining: %d seconds \n',countdown);
     pause(1);
end

% Clear screen and show "Starting..." message briefly
Screen('FillRect', window, grey);
Screen('Flip', window);

fprintf('Countdown completed. Experiment starting now!\n');
fprintf('Press any key to stop the program\n');
fprintf('==========================================\n');

%------------------------------------------
% Initialize timing recording
%------------------------------------------

% Start master timer
tic;
masterStartTime = toc;

% Initialize arrays to store timestamps
stimulusStartTimes = [];
intervalStartTimes = [];
cycleBreakStartTimes = [];
orientationSequence = [];
trialSequence = [];
cycleSequence = [];  % Track which cycle each trial belongs to

%------------------------------------------
% Main experimental loop with error handling
%------------------------------------------

try
    % Set high priority for visual stimulation timing precision
    Priority(topPriorityLevel);
    
    % Initial flip to start timing
    vbl = Screen('Flip', window);
    startTime = vbl;

while ~KbCheck && trialNumber <= maxTrials
    % Get current orientation from the randomized sequence for this cycle
    currentOrientation = currentCycleOrientations(currentOrientationIndex);
    
    %% Visual Stimulation Phase (2 seconds)
    fprintf('  Visual stimulation phase started...\n');
    
    % Record stimulus start time
    stimulusStartTime = toc;
    stimulusStartTimes = [stimulusStartTimes; stimulusStartTime];
    orientationSequence = [orientationSequence; currentOrientation];
    trialSequence = [trialSequence; trialNumber];
    cycleSequence = [cycleSequence; currentCycle];
    
    fprintf('  [%.3f s] Cycle %d, Stimulus %d/%d: Orientation = %d degrees\n', ...
        stimulusStartTime, currentCycle, currentOrientationIndex, numOrientations, currentOrientation);
    
    % Reset phase for each trial
    currentPhase = phase;
    
    % Create properties matrix for current orientation
    propertiesMat = [currentPhase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
    
    % Record start time for this stimulation phase
    stimStartTime = vbl;
    
    % Visual stimulation loop
    frameCount = 0;
    while (vbl - stimStartTime) < stimDuration
        frameCount = frameCount + 1;
        
        % Draw the Gabor with current orientation
        Screen('DrawTextures', window, gabortex, [], [], currentOrientation, [], [], [], [],...
            [], propertiesMat');
        
        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % Update the phase for animation
        propertiesMat(1) = propertiesMat(1) + phasePerFrame;
        
        % Periodically allow other processes (including camera) to update
        % This reduces camera freezing during intensive visual stimulation
        if mod(frameCount, 10) == 0  % Every 10 frames
            drawnow limitrate;  % Allow GUI updates without blocking too long
        end
        
        % Check for key press to exit
        if KbCheck
            break;
        end
    end
    
    % Grey Background Phase (4 seconds)
    fprintf('  Grey background phase started...\n');
    
    % Record interval start time
    intervalStartTime = toc;
    intervalStartTimes = [intervalStartTimes; intervalStartTime];
    
    fprintf('  [%.3f s] Interval started after stimulus %d\n', intervalStartTime, trialNumber);
    
    % Clear screen to grey background
    Screen('FillRect', window, grey);
    vbl = Screen('Flip', window);
    greyStartTime = vbl;
    
    % Grey background loop
    frameCount = 0;
    while (vbl - greyStartTime) < greyDuration
        frameCount = frameCount + 1;
        
        % Keep grey background
        Screen('FillRect', window, grey);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % Periodically allow camera and other GUI updates
        if mod(frameCount, 5) == 0  % More frequent during grey period
            drawnow limitrate;
        end
        
        % Check for key press to exit
        if KbCheck
            break;
        end
    end
    
    % Check for key press to exit (final check)
    if KbCheck
        break;
    end
    
    % Move to next orientation in the randomized sequence
    currentOrientationIndex = currentOrientationIndex + 1;
    
    % Check if we completed a full orientation cycle
    if currentOrientationIndex > numOrientations
        currentOrientationIndex = 1;  % Reset to first orientation
        
        fprintf('  Complete orientation cycle %d/%d finished.\n', currentCycle, maxCycles);
        
        % Increment cycle counter
        currentCycle = currentCycle + 1;
        
        % Check if we've completed all cycles
        if currentCycle > maxCycles
            fprintf('  All %d cycles completed! Experiment will end.\n', maxCycles);
            break; % Exit the main loop
        end
        
        % Generate NEW randomized orientation sequence for next cycle
        currentCycleOrientations = baseOrientations(randperm(numOrientations));
        fprintf('  Cycle %d orientation order: %s\n', currentCycle, mat2str(currentCycleOrientations));
        
        % Add 6-second break between complete orientation cycles
        fprintf('  Starting 6-second break before cycle %d...\n', currentCycle);
        
        % Record cycle break start time
        cycleBreakStartTime = toc;
        cycleBreakStartTimes = [cycleBreakStartTimes; cycleBreakStartTime];
        
        fprintf('  [%.3f s] Cycle break started\n', cycleBreakStartTime);
        
        % Clear screen to grey background for cycle break
        Screen('FillRect', window, grey);
        vbl = Screen('Flip', window);
        cycleBreakStartTime = vbl;
        
        % Cycle break loop
        frameCount = 0;
        while (vbl - cycleBreakStartTime) < cycleBreakDuration
            frameCount = frameCount + 1;
            
            % Keep grey background
            Screen('FillRect', window, grey);
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
            % Allow frequent GUI updates during break
            if mod(frameCount, 3) == 0  % Very frequent during break
                drawnow limitrate;
            end
            
            % Check for key press to exit
            if KbCheck
                break;
            end
        end
        
        fprintf('  Cycle break completed. Starting cycle %d/%d...\n', currentCycle, maxCycles);
    end
    
    % Increment trial number
    trialNumber = trialNumber + 1;
    
    fprintf('  Trial %d completed. Total time: %.2f seconds\n', trialNumber-1, vbl - stimStartTime);
    fprintf('------------------------------------------\n');
    
end

catch ME
    % Handle any errors during the experiment
    fprintf('\nExperiment stopped due to error: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for k = 1:length(ME.stack)
        fprintf('  %s at line %d\n', ME.stack(k).file, ME.stack(k).line);
    end
    
    % Ensure cameras are cleaned up even in case of error
    if exist('cameraInitialized','var') && cameraInitialized
        try
            if exist('backCam','var')
                closePreview(backCam);
                clear backCam;
            end
            if exist('frontCam','var')
                closePreview(frontCam);
                clear frontCam;
            end
            fprintf('Emergency camera cleanup completed.\n');
        catch
            % Ignore cleanup errors
        end
    end
    
    % Reset priority to normal
    try
        Priority(0);
    catch
        % Ignore priority reset errors
    end
    
    % Ensure screen is cleared
    try
        sca;
    catch
        % Ignore screen clearing errors
    end
    
    % Re-throw the error after cleanup
    rethrow(ME);
end

% Calculate total experimental time
totalTime = vbl - startTime;
experimentEndTime = toc;

fprintf('\nExperiment completed!\n');
fprintf('Total cycles completed: %d/%d\n', currentCycle-1, maxCycles);
fprintf('Total trials: %d/%d\n', trialNumber-1, maxTrials);
fprintf('Total time: %.2f seconds (%.2f minutes)\n', totalTime, totalTime/60);
fprintf('Experiment duration (tic/toc): %.3f seconds\n', experimentEndTime);

%------------------------------------------
% Save timing data
%------------------------------------------

% Create timestamp structure
timingData = struct();
timingData.masterStartTime = masterStartTime;
timingData.experimentEndTime = experimentEndTime;
timingData.stimulusStartTimes = stimulusStartTimes;
timingData.intervalStartTimes = intervalStartTimes;
timingData.cycleBreakStartTimes = cycleBreakStartTimes;
timingData.orientationSequence = orientationSequence;
timingData.trialSequence = trialSequence;
timingData.cycleSequence = cycleSequence;
timingData.parameters = struct('stimDuration', stimDuration, ...
                              'greyDuration', greyDuration, ...
                              'cycleBreakDuration', cycleBreakDuration, ...
                              'baseOrientations', baseOrientations, ...
                              'maxCycles', maxCycles, ...
                              'completedCycles', currentCycle-1, ...
                              'maxTrials', maxTrials, ...
                              'randomizedOrder', true);

% Generate filename with timestamp
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('visual_stim_timing_randomized_%s.mat', timeStr);

% Save timing data
save(filename, 'timingData');

fprintf('\n==========================================\n');
fprintf('Timing data saved to: %s\n', filename);
fprintf('==========================================\n');

% Display timing summary
fprintf('\nTiming Summary:\n');
fprintf('Cycles completed: %d/%d\n', currentCycle-1, maxCycles);
fprintf('Number of stimuli presented: %d/%d\n', length(stimulusStartTimes), maxTrials);
fprintf('Number of intervals: %d\n', length(intervalStartTimes));
fprintf('Number of cycle breaks: %d\n', length(cycleBreakStartTimes));

if ~isempty(stimulusStartTimes)
    fprintf('First stimulus started at: %.3f seconds\n', stimulusStartTimes(1));
    fprintf('Last stimulus started at: %.3f seconds\n', stimulusStartTimes(end));
end

% Display orientation presentation count
fprintf('\nOrientation presentation summary:\n');
for ori = baseOrientations
    count = sum(orientationSequence == ori);
    fprintf('  %d degrees: presented %d times\n', ori, count);
end

% Clean up cameras before closing
if cameraInitialized
    try
        % Stop camera previews
        closePreview(backCam);
        closePreview(frontCam);
        % Clear camera objects
        clear backCam frontCam;
        fprintf('Cameras cleaned up successfully.\n');
    catch ME
        warning('Camera cleanup failed: %s', ME.message);
    end
end

% Reset priority to normal before closing
Priority(0);

% Close figure window
if exist('f','var') && isvalid(f)
    close(f);
end

% Clear screen
sca;

% Program completed
fprintf('\nExperiment completed successfully!\n');
fprintf('Data has been saved and screen cleared.\n');
