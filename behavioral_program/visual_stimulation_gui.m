% visual_stimulation_gui.m
%
% A MATLAB program using Psychtoolbox to present drifting gratings
% with parameters set via a graphical user interface (GUI).
% Modified: full-screen sinusoidal grating with smooth sinusoidal fade edges
% and proper start control.

% Clean up workspace and screen
sca;
close all;
clearvars;

% --- Setup Psychtoolbox and Screen ---
PsychDefaultSetup(2);
Screen('Preference','SkipSyncTests',0); % Ensure proper sync in real experiments
h.screenNumber = max(Screen('Screens')); % Use external screen if available

h.white = WhiteIndex(h.screenNumber);
h.grey  = h.white / 2;

% Open a window
[h.window, h.windowRect] = PsychImaging('OpenWindow', h.screenNumber, h.grey);

% --- GUI for parameter input (after screen is initialized) ---
prompt = {
    'Stimulus Duration (s):', ...
    'Inter-Stimulus Interval (ISI) (s):', ...
    'Spatial Frequency (cycles/degree):', ...
    'Orientations (degrees, space-separated):', ...
    'Number of Repeats:', ...
    'Mouse ID:'
    };
dlgtitle = 'Visual Stimulation Parameters';
dims = [1 50];
definput = {
    '4', ...
    '8', ...
    '0.05', ...
    '0 30 60 90 120 150', ...
    '10', ...
    'M1'
    };
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Exit if user cancels
if isempty(answer)
    disp('Program cancelled by user.');
    sca; % Close the PTB window
    return;
end

% Parse GUI input
stimDuration = str2double(answer{1});
isiDuration = str2double(answer{2});
spatialFrequency_cpd = str2double(answer{3});
orientations = str2num(answer{4}); %#ok<ST2NM>
repeats = str2double(answer{5});
mouseID = answer{6};

h.ifi = Screen('GetFlipInterval', h.window);
h.topPriorityLevel = MaxPriority(h.window);
Priority(h.topPriorityLevel);

% Get window dimensions
[h.width, h.height] = Screen('WindowSize', h.window);

% --- Grating Stimulus Setup ---
viewingDistance_cm = 12.5;
monitorWidth_cm = 21; % example
pixels_per_cm = h.width / monitorWidth_cm;
cm_per_degree = 2 * viewingDistance_cm * tan(pi/360);
pixels_per_degree = pixels_per_cm * cm_per_degree;

% Convert spatial frequency from cycles/degree to cycles/pixel
spatialFrequency_cpp = spatialFrequency_cpd / pixels_per_degree;

% Grating parameters
h.gaborDimPix = max(h.width, h.height) + 200; % Make it large enough to cover the screen
h.sigma = min(h.width, h.height) / 6; % Sigma for Gaussian envelope, creating smooth edges
h.contrast = 1.0;
h.aspectRatio = 1.0;
h.phase = 0;
h.phasePerFrame = (360 * 1.5) * h.ifi; % Temporal frequency = 1.5 Hz

% Create procedural gabor texture with proper sinusoidal edge gradient
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;
h.gabortex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

% --- Trial Structure Setup ---
trial_orientations = repmat(orientations, 1, repeats);
trial_sequence = trial_orientations(randperm(length(trial_orientations)));
totalTrials = length(trial_sequence);

% --- Data Saving Setup ---
results.mouseID = mouseID;
results.parameters = struct(...
    'stimDuration', stimDuration, ...
    'isiDuration', isiDuration, ...
    'spatialFrequency_cpd', spatialFrequency_cpd, ...
    'orientations', orientations, ...
    'repeats', repeats ...
);
results.trialLog = cell(totalTrials, 3); % Trial#, Orientation, Timestamp
results.filename = sprintf('visual_stim_log_%s_%s.mat', mouseID, datestr(now, 'yyyymmdd_HHMMSS'));

% --- Start Experiment ---
% Display message in PTB window and wait for keypress
DrawFormattedText(h.window, 'Press any key to start the experiment.', 'center', 'center', h.white);
Screen('Flip', h.window);
KbWait; % Wait for a key press
Screen('Flip', h.window); % Clear the text and show grey screen
WaitSecs(0.5); % Brief pause before starting trials

% Main experiment loop
for trialNum = 1:totalTrials
    currentOrientation = trial_sequence(trialNum);
    fprintf('Trial %d/%d: Orientation = %d degrees\n', trialNum, totalTrials, currentOrientation);

    % Set properties matrix for this trial
    propertiesMat = [h.phase, spatialFrequency_cpp, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];
    
    vbl = Screen('Flip', h.window); 
    startTime = vbl;
    
    while vbl < startTime + stimDuration
        % Draw drifting gabor with sinusoidal edges
        Screen('DrawTextures', h.window, h.gabortex, [], [], currentOrientation, [], [], [], [],...
            kPsychDontDoRotation, propertiesMat');

        % Flip to the screen
        vbl = Screen('Flip', h.window, vbl + 0.5*h.ifi);

        % Update phase for drifting effect
        h.phase = h.phase + h.phasePerFrame;
        propertiesMat(1) = h.phase;
    end
    
    % Log trial data
    results.trialLog{trialNum, 1} = trialNum;
    results.trialLog{trialNum, 2} = currentOrientation;
    results.trialLog{trialNum, 3} = startTime;

    % --- Inter-Stimulus Interval (ISI) ---
    fprintf('ISI period for %f seconds...\n', isiDuration);
    Screen('Flip', h.window); % Show grey screen
    WaitSecs(isiDuration);
    
    % Save results incrementally
    save(results.filename, 'results');
end

% --- End of Experiment ---
DrawFormattedText(h.window, 'Experiment finished!', 'center', 'center', h.white);
Screen('Flip', h.window);
WaitSecs(2);

% Clean up
sca;
Priority(0);
disp('Experiment finished and data saved.');
fprintf('Results saved to: %s\n', results.filename);
