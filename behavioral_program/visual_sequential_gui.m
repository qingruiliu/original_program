%% Visual Sequential GUI - based on stage2 framework
% Modified from stage2_RandITI240116.m 
% Function: 4 seconds visual stimulation with 8 seconds interval
% Random grating orientations presentation

timer = timerfindall;
delete(timer)
sca  
clc
clear h.a
clear all

%% open the monitor with gray background
global h 
PsychDefaultSetup(2);
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

%% GUI for parameter input (similar to visual_stimulation_gui)
prompt = {
    'Stimulus Duration (s):', ...
    'Inter-Stimulus Interval (ISI) (s):', ...
    'Spatial Frequency (cycles/degree):', ...
    'Orientations (degrees, space-separated):', ...
    'Number of Repeats:', ...
    'Mouse ID:'
    };
dlgtitle = 'Visual Sequential Stimulation Parameters';
dims = [1 50];
definput = {
    '4', ...
    '8', ...
    '0.05', ...
    '0 30 60 90 120 150 180 210 240 270 300 330', ...
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
h.stimDuration = str2double(answer{1});
h.isiDuration = str2double(answer{2});
spatialFrequency_cpd = str2double(answer{3});
orientations = str2num(answer{4}); %#ok<ST2NM>
repeats = str2double(answer{5});
h.mouseID = answer{6};

%% Gabor presetting (exactly like stage2)
%size of the gabor patch, full of the height in this case
h.gaborDimPix = h.windowRect(4)*2;
h.width = h.windowRect(3);
h.height = h.windowRect(4);

%center of diplay position
h.center = [(h.width-h.height)/2,0,(h.width+h.height)/2,h.height];

%other parameters
h.sigma = h.gaborDimPix;
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

%% Trial Structure Setup
trial_orientations = repmat(orientations, 1, repeats);
h.trial_sequence = trial_orientations(randperm(length(trial_orientations)));
h.totalTrials = length(h.trial_sequence);

%% program variables (like stage2)
h.data1 = zeros(h.totalTrials,4);

%% Data Saving Setup
results.mouseID = h.mouseID;
results.parameters = struct(...
    'stimDuration', h.stimDuration, ...
    'isiDuration', h.isiDuration, ...
    'spatialFrequency_cpd', spatialFrequency_cpd, ...
    'orientations', orientations, ...
    'repeats', repeats ...
);
results.trialLog = cell(h.totalTrials, 3); % Trial#, Orientation, Timestamp
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results.filename = sprintf('visual_sequential_log_%s_%s.mat', h.mouseID, timestamp);
h.results = results;

%% counterUI (exactly like stage2)
screenSize = get(0,'Screensize'); screenSize(3) = screenSize(3)/2;       
f = figure('Name','Visual Sequential Monitor','Position',screenSize,'Color',[0.95 0.95 0.95]);          %open the monitor UI

%total counter UI
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.05 0.95 0.1 0.04],'String','Trial Number','BackgroundColor',[0 1 1],...
    'FontSize',14);
h.totalTrialNumUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.05 0.91 0.1 0.03],'FontSize',14); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.18 0.95 0.12 0.04],'String','Current Orientation','BackgroundColor',[0 1 1],...
    'FontSize',14);
h.orientationUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.18 0.91 0.12 0.03],'FontSize',14); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.32 0.95 0.1 0.04],'String','Progress','BackgroundColor',[1 1 0],...
    'FontSize',14);
h.progressUI= uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.32 0.91 0.1 0.03],'FontSize',14); 

%visual stimulation parameters
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.11 0.86 0.25 0.03],'String','Stimulation Parameters','BackgroundColor',[0.83 0.5 1],...
    'FontSize',14);

VSPara = {'Stim Duration' 'ISI Duration' 'Total Trials'};
VSParaValue = [h.stimDuration h.isiDuration h.totalTrials];

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.05 0.81 0.12 0.04],'String',VSPara{1},'BackgroundColor',[1 1 1],...
    'FontSize',14);
h.stimDurUI =uicontrol(f,'Style','edit','String',num2str(VSParaValue(1)),'Units','normalized',...
    'Position',[0.05 0.79 0.12 0.03],'FontSize',14); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.18 0.81 0.12 0.04],'String',VSPara{2},'BackgroundColor',[1 1 1],...
    'FontSize',14);
h.isiDurUI =uicontrol(f,'Style','edit','String',num2str(VSParaValue(2)),'Units','normalized',...
    'Position',[0.18 0.79 0.12 0.03],'FontSize',14); 

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.31 0.81 0.12 0.04],'String',VSPara{3},'BackgroundColor',[1 1 1],...
    'FontSize',14);
h.totalTrialsUI =uicontrol(f,'Style','edit','String',num2str(VSParaValue(3)),'Units','normalized',...
    'Position',[0.31 0.79 0.12 0.03],'FontSize',14); 

% Status display
uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.5 0.95 0.15 0.04],'String','Current Status','BackgroundColor',[1 1 1],...
    'FontSize',14);
h.statusUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.5 0.9 0.15 0.04],'FontSize',14);

uicontrol(f,'Style','text','Units','normalized',...
    'Position',[0.7 0.95 0.15 0.04],'String','Elapsed Time','BackgroundColor',[1 1 1],...
    'FontSize',14);
h.elapsedTimeUI = uicontrol(f,'Style','edit','Units','normalized',...
    'Position',[0.7 0.9 0.15 0.04],'FontSize',14);

% trial raster (like stage2)
h.trialRaster = axes(f,'Position',[0.1 0.44 0.3 0.3],'FontSize',14);
title('Trial Timeline');
xlabel('Time (seconds)');
ylabel('Trial Number');
xlim([0 h.stimDuration + h.isiDuration]);
ylim([0 h.totalTrials]);
set(h.trialRaster,'Ydir','reverse');
hold on

% time length monitor (like stage2)
h.trialTime = axes(f,'Position',[0.55 0.44 0.3 0.3],'FontSize',14);
title('Trial Duration');
xlabel('Trial Number');     
xlim([1 h.totalTrials]);
ylabel('Duration (seconds)');
hold on

%open the cameras (exactly like stage2)
% Get camera resolutions
h.backCam = webcam(1);
h.frontCam = webcam(2);

backRes = str2double(strsplit(h.backCam.Resolution,'x'));
frontRes = str2double(strsplit(h.frontCam.Resolution,'x'));

% Calculate aspect ratios
backAspect = backRes(1) / backRes(2);
frontAspect = frontRes(1) / frontRes(2);

% Set UI positions based on aspect ratios (normalized units)
backHeight = 0.25;
backWidth = backAspect * backHeight;
frontHeight = 0.25;
frontWidth = frontAspect * frontHeight;

% Place back camera UI
h.backCamUI = axes(f, 'Position', [0.55 0.1 backWidth backHeight]);
h.im = image(zeros(backRes(2), backRes(1), 3, 'uint8'), 'Parent', h.backCamUI);
preview(h.backCam, h.im);
text(h.backCamUI, 0.5, -0.1, 'Back Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 14);  % Add label below the image

% Place front camera UI
h.frontCamUI = axes(f, 'Position', [0.05 0.1 frontWidth frontHeight]);
title(h.frontCamUI, 'Front Camera');
h.im2 = image(zeros(frontRes(2), frontRes(1), 3, 'uint8'), 'Parent', h.frontCamUI);
preview(h.frontCam, h.im2);
text(h.frontCamUI, 0.5, -0.1, 'Front Camera', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 14);  % Add label below the image

%% 1 minute countdown (modified from stage2)
countDown;

%% start tic (like stage2)
totalTic = tic;
h.trialNum = [];
experimentStartTime = GetSecs;

%% main loop (modified from stage2)
for trialNum = 1:h.totalTrials
    h.trialNum = trialNum;
    h.currentOrientation = h.trial_sequence(trialNum);
    
    % Process display rotation (same as visual_stimulation_gui.m)
    rotatedOrientation = mod(h.currentOrientation + 180, 360);
    if rotatedOrientation == 0 || rotatedOrientation == 360
        h.displayOrientation = 0;
    else
        h.displayOrientation = 360 - rotatedOrientation;
    end
    
    h.trialGlobalTic = tic;
    disp('------------------new trial-----------------------')
    fprintf('trialNum = %s, Orientation = %d°\n',num2str(trialNum), h.currentOrientation)
    set(h.totalTrialNumUI,'String',num2str(trialNum));
    set(h.orientationUI,'String',sprintf('%d°', h.currentOrientation));
    progressPercent = (trialNum-1) / h.totalTrials * 100;
    set(h.progressUI,'String',sprintf('%.1f%%', progressPercent));

    visiStim;      %present the visual stimulation
    
    ISIperiod;      %inter-stimulus interval

    fprintf('Trial %d completed: Orientation = %d°\n', trialNum, h.currentOrientation);
    
    % Record trial data
    h.results.trialLog{trialNum, 1} = trialNum;
    h.results.trialLog{trialNum, 2} = h.currentOrientation;
    h.results.trialLog{trialNum, 3} = round(GetSecs - experimentStartTime, 4);
    
    h.data1(trialNum,1) = trialNum;
    h.data1(trialNum,2) = h.currentOrientation;
    h.data1(trialNum,3) = toc(h.trialGlobalTic); %time length for individual trials
    h.data1(trialNum,4) = GetSecs - experimentStartTime;
    
    trialTimeLength = h.data1(trialNum,3);
    plot(h.trialTime,h.data1(trialNum,1),h.data1(trialNum,3),'-ok');
    
    % Plot trial timeline
    plot(h.trialRaster, [0 h.stimDuration], [trialNum trialNum], 'g-', 'LineWidth', 3); % Stimulus period
    plot(h.trialRaster, [h.stimDuration h.stimDuration+h.isiDuration], [trialNum trialNum], 'b-', 'LineWidth', 1); % ISI period
    
    % Update elapsed time
    elapsedTime = toc(totalTic);
    set(h.elapsedTimeUI,'String',sprintf('%.1f min', elapsedTime/60));
    
    % Save data periodically (every 20 trials)
    if mod(trialNum, 20) == 0
        try
            save(h.results.filename, 'h');
            fprintf('Data saved at trial %d\n', trialNum);
        catch ME
            warning('Failed to save data at trial %d: %s', trialNum, ME.message);
        end
    end
end

%% End of experiment (like stage2)
set(h.statusUI,'String','Experiment Finished');
set(h.progressUI,'String','100%');

disp('---------------------finished: all trials!---------------')
totalTime = toc(totalTic);

% Final data save
try
    save(h.results.filename, 'h');
    fprintf('Final data save completed successfully\n');
catch ME
    warning('Final data save failed: %s', ME.message);
    backup_filename = sprintf('backup_%s', h.results.filename);
    try
        save(backup_filename, 'h');
        fprintf('Data saved to backup file: %s\n', backup_filename);
    catch
        warning('Backup save also failed. Data may be lost.');
    end
end

% Cleanup cameras (like stage2)
stoppreview(h.backCam);
clear h.backCam;
stoppreview(h.frontCam);
clear h.frontCam;

sca;
fprintf('>> total time cost:  %s minutes %s seconds \n',num2str(floor((totalTime)/60)),num2str(mod(totalTime,60)));
fprintf('Results saved to: %s\n', h.results.filename);
h.data1
save h

%% functions (modified from stage2)

function updateVbl(~,~)
    global h
    h.vbl = Screen('Flip',h.window);
    h.vblt0 = h.vbl;
end

function countDown(~,~)
    countdownDuration = 10;
    disp('---------Countdown started...check the setup!!!--------------------')

    for remainingSeconds = countdownDuration :-1 :0
        fprintf('Time remaining: %d seconds \n',remainingSeconds);
        pause(1);
    end
    disp('Countdown finished! Start visual stimulation!')
end

%% trial procedure functions (modified from stage2)

function visiStim(~,~)
    global h
    set(h.statusUI,'String','Visual Stimulation');
    updateVbl;
    
    % Reset phase (similar to stage2)
    h.propertiesMat(1) = 0;
    
    % Stimulus presentation loop (based on stage2 visiStim)
    while h.vbl - h.vblt0 <= h.stimDuration
        Screen('DrawTextures', h.window, h.gabortex, [], [], h.displayOrientation, [], [], [], [],...
            kPsychDontDoRotation, h.propertiesMat');
        h.vbl = Screen('Flip', h.window, h.vbl + (h.waitframes - 0.5) * h.ifi);   %double buffering 
        h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
    end
    updateVbl;
    disp('>>Visual stimulation ended!! ISI starts!!')
end

function ISIperiod(~,~)
    global h
    set(h.statusUI,'String','Inter-Stimulus Interval');
    h.ISITime = tic;
    while toc(h.ISITime) <= h.isiDuration    
    end
    fprintf('>> ISI period finished!!! Time length is %s \n',num2str(h.isiDuration))
end
