% Extended Visual Stimulation Program with Random TF and SF
% Displays 2 seconds of visual stimulation followed by 4 seconds of grey background
% Orientation changes by 30 degrees each cycle
% Random Temporal Frequency (TF) and Spatial Frequency (SF) for each trial
% Clear the workspace and the screen
clc;
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

%--------------------
% Physical setup parameters
%--------------------
% Monitor specifications
monitorWidth_mm = 210.4;   % Visible area width in mm
monitorHeight_mm = 157.8;  % Visible area height in mm
resolutionX = 1024;        % Horizontal resolution in pixels
resolutionY = 768;         % Vertical resolution in pixels
refreshRate = 60;          % Monitor refresh rate in Hz
viewingDistance_cm = 20;   % Distance from mouse eye to screen in cm

% Calculate visual angle per pixel
viewingDistance_mm = viewingDistance_cm * 10;  % Convert to mm
pixelSize_mm = monitorWidth_mm / resolutionX;  % Size of one pixel in mm
visualAngle_perPixel_rad = 2 * atan(pixelSize_mm / (2 * viewingDistance_mm));  % in radians
visualAngle_perPixel_deg = rad2deg(visualAngle_perPixel_rad);  % in degrees

fprintf('\n========== Visual Angle Calculation ==========\n');
fprintf('Monitor size: %.1f × %.1f mm\n', monitorWidth_mm, monitorHeight_mm);
fprintf('Resolution: %d × %d pixels\n', resolutionX, resolutionY);
fprintf('Viewing distance: %.1f cm\n', viewingDistance_cm);
fprintf('Pixel size: %.4f mm\n', pixelSize_mm);
fprintf('Visual angle per pixel: %.4f degrees\n', visualAngle_perPixel_deg);
fprintf('==============================================\n\n');

% Temporal Frequency (TF) values in Hz
temporalFrequencies = [0.5, 1, 2, 4, 8];  % 5 different TF values
numTF = length(temporalFrequencies);

% Spatial Frequency (SF) values in cycles per degree (cpd)
% Common SF values for mouse vision research
spatialFrequencies_cpd = [0.02, 0.04, 0.08, 0.016, 0.32];  % cycles per degree
numSF = length(spatialFrequencies_cpd);

% Convert SF from cycles/degree to cycles/pixel
spatialFrequencies_cpp = spatialFrequencies_cpd * visualAngle_perPixel_deg;  % cycles per pixel

fprintf('Spatial Frequencies:\n');
for i = 1:numSF
    fprintf('  SF %d: %.4f cpd = %.6f cycles/pixel\n', i, spatialFrequencies_cpd(i), spatialFrequencies_cpp(i));
end
fprintf('\n');

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

% Orientation parameters - 30 degree increments (12 total orientations)
orientations = 0:30:330;  % [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330] degrees
numOrientations = length(orientations);
currentOrientationIndex = 1;

% Session and trial parameters
numSessions = 5;  % 总共5个session
trialsPerSession = 60;  % 每个session有60个trial
maxTrials = numSessions * trialsPerSession;  % 总试次数 = 5 × 60 = 300
currentSession = 1;  % 当前session数
trialInSession = 1;  % 当前session中的trial数

% Cycle and trial parameters
maxCycles = 25;  % 总共25个循环 (25 × 12 = 300 trials)
currentCycle = 1;  % 当前循环数

% Trial counter
trialNumber = 1;

% Generate random TF and SF sequence for all trials
% Each trial gets a random combination of TF and SF
rng('shuffle');  % Initialize random number generator with current time
randomTFIndex = randi(numTF, maxTrials, 1);
randomSFIndex = randi(numSF, maxTrials, 1);
randomTF = temporalFrequencies(randomTFIndex);
randomSF_cpd = spatialFrequencies_cpd(randomSFIndex);  % In cycles/degree
randomSF_cpp = spatialFrequencies_cpp(randomSFIndex);  % In cycles/pixel

% Display initial information
fprintf('Extended Visual Stimulation Program with Random TF and SF\n');
fprintf('Orientations: %s degrees (30-degree increments)\n', mat2str(orientations));
fprintf('Temporal Frequencies: %s Hz\n', mat2str(temporalFrequencies));
fprintf('Spatial Frequencies: %s cycles/degree\n', mat2str(spatialFrequencies_cpd));
fprintf('Total sessions: %d\n', numSessions);
fprintf('Trials per session: %d\n', trialsPerSession);
fprintf('Total cycles: %d (12 stimuli per cycle)\n', maxCycles);
fprintf('Total trials: %d\n', maxTrials);
fprintf('Stimulation duration: %d seconds\n', stimDuration);
fprintf('Grey background duration: %d seconds\n', greyDuration);
fprintf('Cycle break duration: %d seconds\n', cycleBreakDuration);
fprintf('==========================================\n');
%% create GUI
screenSize = get(0,'Screensize'); 
screenSize(3) = screenSize(3)/2;         
f = figure('Name','trial monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);   

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

    % Set UI positions based on aspect ratios (normalized units)
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
    title(frontCamUI, 'Back Camera');
    im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', frontCamUI);
    preview(frontCam, im2);
    text(frontCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    
    cameraInitialized = true;
    fprintf('Camera initialization successful.\n');
    
    % Lower priority slightly to allow camera operation
    Priority(0);  % Reset to normal priority after camera setup
catch ME
    warning('Camera initialization failed: %s. Continuing without cameras.', ME.message);
    cameraInitialized = false;
end

%% Add monitoring display to GUI
% Create monitoring panel
monitorPanel = uipanel(f, 'Title', 'Experiment Monitor', ...
    'Position', [0.05 0.65 0.35 0.28], ...
    'BackgroundColor', [0.95 0.95 0.95], ...
    'FontSize', 14, 'FontWeight', 'bold');

% Session display
sessionText = uicontrol(monitorPanel, 'Style', 'text', ...
    'String', sprintf('Current Session: %d / %d', currentSession, numSessions), ...
    'Units', 'normalized', 'Position', [0.1 0.7 0.8 0.2], ...
    'FontSize', 16, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95 0.95 0.95], ...
    'HorizontalAlignment', 'left');

% Trial in session display
trialText = uicontrol(monitorPanel, 'Style', 'text', ...
    'String', sprintf('Trial in Session: %d / %d', trialInSession, trialsPerSession), ...
    'Units', 'normalized', 'Position', [0.1 0.45 0.8 0.2], ...
    'FontSize', 16, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95 0.95 0.95], ...
    'HorizontalAlignment', 'left');

% Status display
statusText = uicontrol(monitorPanel, 'Style', 'text', ...
    'String', 'Status: Waiting to start', ...
    'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.2], ...
    'FontSize', 16, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95 0.95 0.95], ...
    'ForegroundColor', [0 0.5 0], ...
    'HorizontalAlignment', 'left');

drawnow;
%%
%------------------------------------------
% Confirmation and countdown before experiment
%------------------------------------------
msg = msgbox('Ready to start Session 1?');
waitfor(msg)
fprintf('Starting 30-second countdown for Session 1...\n');
set(statusText, 'String', 'Status: Countdown to Session 1', 'ForegroundColor', [0 0.5 0]);
drawnow;

% 30-second countdown
for countdown = 30:-1:1
    fprintf('Time remaining: %d seconds \n', countdown);
    set(statusText, 'String', sprintf('Status: Starting in %d seconds', countdown), ...
        'ForegroundColor', [0 0.5 0]);
    drawnow;
    pause(1);
end

% Clear screen and show "Starting..." message briefly
Screen('FillRect', window, grey);
Screen('Flip', window);

set(statusText, 'String', 'Status: Running', 'ForegroundColor', [0 0.7 0]);
drawnow;

fprintf('Countdown completed. Experiment starting now!\n');
fprintf('Press any key to stop the program\n');
fprintf('==========================================\n');

%------------------------------------------
% Initialize timing recording
%------------------------------------------

% Start master timer
tic;
masterStartTime = toc;

% Initialize arrays to store timestamps and parameters
stimulusStartTimes = [];
intervalStartTimes = [];
cycleBreakStartTimes = [];
sessionBreakStartTimes = [];  % Store session break start times
orientationSequence = [];
trialSequence = [];
sessionSequence = [];  % Store which session each trial belongs to
tfSequence = [];
sfSequence_cpd = [];
sfSequence_cpp = [];

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
    % Update GUI monitoring display
    set(sessionText, 'String', sprintf('Current Session: %d / %d', currentSession, numSessions));
    set(trialText, 'String', sprintf('Trial in Session: %d / %d', trialInSession, trialsPerSession));
    set(statusText, 'String', 'Status: Running', 'ForegroundColor', [0 0.7 0]);
    drawnow limitrate;
    
    % Get current orientation
    currentOrientation = orientations(currentOrientationIndex);
    
    % Get random TF and SF for this trial
    currentTF = randomTF(trialNumber);       % Temporal frequency in Hz
    currentSF_cpd = randomSF_cpd(trialNumber);  % Spatial frequency in cycles/degree
    currentSF_cpp = randomSF_cpp(trialNumber);  % Spatial frequency in cycles/pixel
    
    % Calculate phase increment per frame based on temporal frequency
    % For PTB CreateProceduralSineGrating: phase is in DEGREES (0-360)
    % TF (Hz) = cycles per second = 360 degrees per second
    % phasePerFrame = (TF * 360 degrees) / frameRate (degrees per frame)
    frameRate = 1 / ifi;
    phasePerFrame = (currentTF * 360) / frameRate;  % degrees per frame
    
    %% Visual Stimulation Phase (2 seconds)
    fprintf('  Visual stimulation phase started...\n');
    
    % Record stimulus start time and parameters
    stimulusStartTime = toc;
    stimulusStartTimes = [stimulusStartTimes; stimulusStartTime];
    orientationSequence = [orientationSequence; currentOrientation];
    trialSequence = [trialSequence; trialNumber];
    sessionSequence = [sessionSequence; currentSession];  % Record session number
    tfSequence = [tfSequence; currentTF];
    sfSequence_cpd = [sfSequence_cpd; currentSF_cpd];
    sfSequence_cpp = [sfSequence_cpp; currentSF_cpp];
    
    fprintf('  [%.3f s] Trial %d: Ori = %d°, TF = %.1f Hz, SF = %.4f cpd (%.6f cyc/pix)\n', ...
        stimulusStartTime, trialNumber, currentOrientation, currentTF, currentSF_cpd, currentSF_cpp);

    % Reset phase for each trial
    currentPhase = phase;
    
    % Create properties matrix for current orientation and spatial frequency
    % Use cycles/pixel for the Gabor texture
    propertiesMat = [currentPhase, currentSF_cpp, sigma, contrast, aspectRatio, 0, 0, 0];
    
    % Record start time for this stimulation phase
    stimStartTime = vbl;
    
    % Visual stimulation loop
    frameCount = 0;
    while (vbl - stimStartTime) < stimDuration
        frameCount = frameCount + 1;
        
        % Draw the Gabor with current orientation and parameters
        Screen('DrawTextures', window, gabortex, [], [], currentOrientation, [], [], [], [],...
            [], propertiesMat');
        
        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        % Update the phase for animation based on temporal frequency
        propertiesMat(1) = propertiesMat(1) + phasePerFrame;
        
        % Periodically allow other processes (including camera) to update
        if mod(frameCount, 10) == 0  % Every 10 frames
            drawnow limitrate;
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
    
    fprintf('  [%.3f s] Interval started after trial %d\n', intervalStartTime, trialNumber);
    
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
    
    % Move to next orientation
    currentOrientationIndex = currentOrientationIndex + 1;
    
    % Check if we completed a full orientation cycle
    if currentOrientationIndex > numOrientations
        currentOrientationIndex = 1;  % Reset to first orientation
        
        % Increment cycle counter
        currentCycle = currentCycle + 1;
        
        fprintf('  Complete orientation cycle %d/%d finished.\n', currentCycle-1, maxCycles);
        
        % Check if we've completed all cycles
        if currentCycle > maxCycles
            fprintf('  All %d cycles completed! Experiment will end.\n', maxCycles);
            break; % Exit the main loop
        end
        
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
    trialInSession = trialInSession + 1;
    
    fprintf('  Trial %d (Session %d, Trial %d) completed. Total time: %.2f seconds\n', ...
        trialNumber-1, currentSession, trialInSession-1, vbl - stimStartTime);
    fprintf('------------------------------------------\n');
    
    % Check if session is complete
    if trialInSession > trialsPerSession && currentSession < numSessions
        fprintf('\n========== Session %d Complete ==========\n', currentSession);
        set(statusText, 'String', sprintf('Status: Session %d Complete!', currentSession), ...
            'ForegroundColor', [0 0 1]);
        drawnow;
        
        % Record session break start time
        sessionBreakStartTime = toc;
        sessionBreakStartTimes = [sessionBreakStartTimes; sessionBreakStartTime];
        fprintf('[%.3f s] Session break started\n', sessionBreakStartTime);
        
        % Show confirmation dialog
        msg = msgbox(sprintf('Session %d completed! Ready for Session %d?', ...
            currentSession, currentSession + 1));
        waitfor(msg);
        
        % 30-second countdown before next session
        fprintf('Starting 30-second countdown before Session %d...\n', currentSession + 1);
        for countdown = 30:-1:1
            fprintf('Time remaining: %d seconds \n', countdown);
            set(statusText, 'String', sprintf('Status: Next session in %d seconds', countdown), ...
                'ForegroundColor', [1 0.5 0]);
            drawnow;
            pause(1);
        end
        
        % Move to next session
        currentSession = currentSession + 1;
        trialInSession = 1;
        
        fprintf('Starting Session %d...\n', currentSession);
        set(statusText, 'String', 'Status: Running', 'ForegroundColor', [0 0.7 0]);
        set(sessionText, 'String', sprintf('Current Session: %d / %d', currentSession, numSessions));
        set(trialText, 'String', sprintf('Trial in Session: %d / %d', trialInSession, trialsPerSession));
        drawnow;
        
        % Clear screen for next session
        Screen('FillRect', window, grey);
        vbl = Screen('Flip', window);
    end
    
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
timingData.sessionBreakStartTimes = sessionBreakStartTimes;  % Session break times
timingData.orientationSequence = orientationSequence;
timingData.trialSequence = trialSequence;
timingData.sessionSequence = sessionSequence;  % Which session each trial belongs to
timingData.tfSequence = tfSequence;  % Store TF sequence (Hz)
timingData.sfSequence_cpd = sfSequence_cpd;  % Store SF sequence (cycles/degree)
timingData.sfSequence_cpp = sfSequence_cpp;  % Store SF sequence (cycles/pixel)
timingData.parameters = struct('stimDuration', stimDuration, ...
                              'greyDuration', greyDuration, ...
                              'cycleBreakDuration', cycleBreakDuration, ...
                              'orientations', orientations, ...
                              'temporalFrequencies', temporalFrequencies, ...
                              'spatialFrequencies_cpd', spatialFrequencies_cpd, ...
                              'spatialFrequencies_cpp', spatialFrequencies_cpp, ...
                              'viewingDistance_cm', viewingDistance_cm, ...
                              'monitorWidth_mm', monitorWidth_mm, ...
                              'monitorHeight_mm', monitorHeight_mm, ...
                              'resolution', [resolutionX, resolutionY], ...
                              'refreshRate', refreshRate, ...
                              'visualAngle_perPixel_deg', visualAngle_perPixel_deg, ...
                              'maxCycles', maxCycles, ...
                              'completedCycles', currentCycle-1, ...
                              'maxTrials', maxTrials, ...
                              'numSessions', numSessions, ...
                              'trialsPerSession', trialsPerSession, ...
                              'completedSessions', currentSession);

% Generate filename with timestamp
timeStr = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('visual_stim_randomTF_SF_timing_%s.mat', timeStr);

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

% Display TF and SF usage statistics
fprintf('\nTemporal Frequency (TF) usage:\n');
for i = 1:numTF
    count = sum(tfSequence == temporalFrequencies(i));
    fprintf('  TF = %.1f Hz: %d trials (%.1f%%)\n', ...
        temporalFrequencies(i), count, 100*count/length(tfSequence));
end

fprintf('\nSpatial Frequency (SF) usage:\n');
for i = 1:numSF
    count = sum(sfSequence_cpd == spatialFrequencies_cpd(i));
    fprintf('  SF = %.4f cpd (%.6f cyc/pix): %d trials (%.1f%%)\n', ...
        spatialFrequencies_cpd(i), spatialFrequencies_cpp(i), count, 100*count/length(sfSequence_cpd));
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
