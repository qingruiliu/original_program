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
Screen('Preference','ScreenToHead',0,0,1);
Screen('Preference','ScreenToHead',1,0,2);
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
% 使用与stage2相同的简单设置方式
h.gaborDimPix = max(h.width, h.height); % 使用屏幕尺寸
h.contrast = 1.0;
h.phase = 0;

% 使用与stage2相同的空间频率设置
h.numCycles = 7; % 与stage2相同
h.freq = h.numCycles / h.gaborDimPix; % 与stage2相同的计算方式

% 使用与stage2相同的参数
h.sigma = h.gaborDimPix;
h.aspectRatio = 1;
h.backgroundOffset = [0.5 0.5 0.5 0.0];
h.disableNorm = 1;
h.preContrastMultiplier = 0.5;

% 创建Gabor纹理而不是正弦光栅，与stage2保持一致
h.gratingtex = CreateProceduralGabor(h.window, h.gaborDimPix, h.gaborDimPix, [],...
    h.backgroundOffset, h.disableNorm, h.preContrastMultiplier);

% 属性矩阵与stage2相同
h.propertiesMat = [h.phase, h.freq, h.sigma, h.contrast, h.aspectRatio, 0, 0, 0];

% 使用与stage2相同的相位增量
h.phasePerFrame = 5 * pi; % 与stage2完全相同
h.waitframes = 1; % 与stage2相同

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
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results.filename = sprintf('visual_stim_log_%s_%s.mat', mouseID, timestamp);

% --- 初始化摄像头 ---
try
    % 初始化刺激监视摄像头
    h.cam1 = webcam(1);
    h.cam1Available = true;
catch
    warning('Stimulus monitoring camera (webcam 1) not available');
    h.cam1Available = false;
end

try
    % 初始化小鼠状态监视摄像头
    h.cam2 = webcam(2);
    h.cam2Available = true;
catch
    warning('Mouse monitoring camera (webcam 2) not available');
    h.cam2Available = false;
end

% --- 创建监视状态GUI ---
h.monitorFig = figure('Name', 'Experiment Monitor', 'Position', [100, 100, 800, 600], ...
    'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');

% 创建两个摄像头视图面板
h.cam1Panel = uipanel('Parent', h.monitorFig, 'Title', 'Back Camera - Stimulus Monitor', ...
    'Position', [0.025, 0.525, 0.45, 0.45], 'FontSize', 12, 'FontWeight', 'bold');

h.cam2Panel = uipanel('Parent', h.monitorFig, 'Title', 'Front Camera - Mouse Monitor', ...
    'Position', [0.525, 0.525, 0.45, 0.45], 'FontSize', 12, 'FontWeight', 'bold');

% 创建摄像头轴
h.cam1Axes = axes('Parent', h.cam1Panel, 'Position', [0.05, 0.05, 0.9, 0.9]);
if h.cam1Available
    h.cam1Image = image(h.cam1Axes, snapshot(h.cam1));
    axis(h.cam1Axes, 'off');
else
    text(h.cam1Axes, 0.5, 0.5, 'Camera 1 not connected', 'HorizontalAlignment', 'center');
    axis(h.cam1Axes, [0 1 0 1]);
    axis(h.cam1Axes, 'off');
end

h.cam2Axes = axes('Parent', h.cam2Panel, 'Position', [0.05, 0.05, 0.9, 0.9]);
if h.cam2Available
    h.cam2Image = image(h.cam2Axes, snapshot(h.cam2));
    axis(h.cam2Axes, 'off');
else
    text(h.cam2Axes, 0.5, 0.5, 'Camera 2 not connected', 'HorizontalAlignment', 'center');
    axis(h.cam2Axes, [0 1 0 1]);
    axis(h.cam2Axes, 'off');
end

% 创建状态显示面板
h.statusPanel = uipanel('Parent', h.monitorFig, 'Title', 'Experiment Status', ...
    'Position', [0.025, 0.05, 0.45, 0.425], 'FontSize', 12, 'FontWeight', 'bold');

h.statusAxes = axes('Parent', h.statusPanel, 'Position', [0.05, 0.05, 0.9, 0.9]);
hold(h.statusAxes, 'on');
axis(h.statusAxes, [0 1 0 1]);
axis(h.statusAxes, 'off');

% 创建方向信息文本
h.orientText = text(h.statusAxes, 0.1, 0.8, 'Direction: --°', 'FontSize', 14);

% 创建试验计数文本
h.trialText = text(h.statusAxes, 0.1, 0.6, 'Trial: 0 / 0', 'FontSize', 14);

% 创建刺激/间隔状态指示器
h.statusText = text(h.statusAxes, 0.1, 0.4, 'Status: Waiting to start', 'FontSize', 14);

% 创建时间指示器
h.timeText = text(h.statusAxes, 0.1, 0.2, 'Time: 0 s', 'FontSize', 14);

% 创建进度条面板
h.progressPanel = uipanel('Parent', h.monitorFig, 'Title', 'Overall Progress', ...
    'Position', [0.525, 0.05, 0.45, 0.425], 'FontSize', 12, 'FontWeight', 'bold');

h.progressBar = axes('Parent', h.progressPanel, 'Position', [0.1, 0.4, 0.8, 0.2]);
h.progressPatch = fill(h.progressBar, [0 0 0 0], [0 0 1 1], 'g');
xlim(h.progressBar, [0 1]);
ylim(h.progressBar, [0 1]);
set(h.progressBar, 'YTick', [], 'XTick', [0 0.25 0.5 0.75 1], ...
    'XTickLabel', {'0%', '25%', '50%', '75%', '100%'});

% 创建更新监视器的函数句柄
h.updateMonitor = @(trialNum, orientation, state, timeElapsed, totalProgress) updateMonitorGUI(h, trialNum, totalTrials, orientation, state, timeElapsed, totalProgress);

% 预先更新一次GUI
h.updateMonitor(0, 0, 'READY', 0, 0);
drawnow;

% --- Start Experiment ---
uiwait(msgbox('Press OK to start the experiment.')); % <-- waits for user

% 添加10秒倒计时功能
fprintf('\n========== countdown start: 10 seconds ========== \n');
for countdown = 10:-1:1
    fprintf('countdown: %d seconds...\n', countdown);
    pause(1);  % 等待1秒
end
fprintf('================================\n\n');

disp('Experiment will start.');

% 记录实验开始时间作为基准时间
experimentStartTime = GetSecs;

% Main experiment loop
for trialNum = 1:totalTrials
    currentOrientation = trial_sequence(trialNum);
    
    % 步骤1: 处理显示器顺时针旋转180度
    rotatedOrientation = mod(currentOrientation + 180, 360);
    
    % 步骤2: 转换从顺时针系统到逆时针系统
    if rotatedOrientation == 0 || rotatedOrientation == 360
        displayOrientation = 0;
    else
        displayOrientation = 360 - rotatedOrientation;
    end
    
    fprintf('Trial %d/%d: Orientation = %d°\n', trialNum, totalTrials, currentOrientation);
    
    % 更新监视器状态为准备
    totalProgress = (trialNum - 1) / totalTrials;
    h.updateMonitor(trialNum, currentOrientation, 'READY', 0, totalProgress);
    
    % 使用与stage2相同的VBL更新方式
    vbl = Screen('Flip', h.window);
    vblt0 = vbl;
    startTime = vbl;
    
    % 重置相位（与stage2类似）
    h.propertiesMat(1) = 0;
    
    % 减少GUI更新频率的计数器
    guiUpdateCounter = 0;
    
    % 刺激呈现循环 - 使用与stage2相同的方式
    while vbl - vblt0 <= stimDuration
        % 每10帧更新一次GUI以减少卡顿
        guiUpdateCounter = guiUpdateCounter + 1;
        if mod(guiUpdateCounter, 10) == 1
            timeElapsed = round(vbl - vblt0);
            h.updateMonitor(trialNum, currentOrientation, 'STIMULUS', timeElapsed, totalProgress + (timeElapsed/stimDuration)/totalTrials/2);
        end
        
        % 使用与stage2相同的绘制方式
        Screen('DrawTexture', h.window, h.gratingtex, [], [], displayOrientation, [], [], [], [],...
            kPsychDontDoRotation, h.propertiesMat');
        
        % 使用与stage2相同的翻转时序
        vbl = Screen('Flip', h.window, vbl + (h.waitframes - 0.5) * h.ifi);
        
        % 使用与stage2相同的相位更新
        h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
    end
    
    % Log trial data - 保存相对于实验开始时间的时间戳（秒）
    results.trialLog{trialNum, 1} = trialNum;
    results.trialLog{trialNum, 2} = currentOrientation;
    results.trialLog{trialNum, 3} = round(startTime - experimentStartTime, 4); % 保留4位小数的相对时间

    % --- Inter-Stimulus Interval (ISI) ---
    fprintf('ISI period for %f seconds...\n', isiDuration);
    
    % 间隔开始时间
    isiStartTime = Screen('Flip', h.window); % Show grey screen
    
    % 间隔期间减少更新频率
    lastUpdateTime = 0;
    while GetSecs < isiStartTime + isiDuration
        currentTime = GetSecs - isiStartTime;
        % 每秒更新一次GUI
        if floor(currentTime) > lastUpdateTime
            timeElapsed = floor(currentTime);
            currentProgress = totalProgress + 0.5/totalTrials + (timeElapsed/isiDuration)/totalTrials/2;
            h.updateMonitor(trialNum, currentOrientation, 'INTERVAL', timeElapsed, currentProgress);
            lastUpdateTime = floor(currentTime);
        end
        WaitSecs(0.1);  % 使用WaitSecs而不是pause以减少系统负载
    end
    
    % Save results incrementally
    save(results.filename, 'results');
end

% --- End of Experiment ---
h.updateMonitor(totalTrials, 0, '实验完成', 0, 1);

%DrawFormattedText(h.window, 'Experiment finished!', 'center', 'center', h.white);
Screen('Flip', h.window);
WaitSecs(2);

% 清理摄像头资源
if h.cam1Available
    clear h.cam1;
end
if h.cam2Available
    clear h.cam2;
end

% Clean up
sca;
Priority(0);
disp('Experiment finished and data saved.');
fprintf('Results saved to: %s\n', results.filename);

function updateMonitorGUI(h, trialNum, totalTrials, orientation, state, timeElapsed, totalProgress)
% updateMonitorGUI - Update the experiment monitoring GUI status
%
%   Parameters:
%     h - Handle structure
%     trialNum - Current trial number
%     totalTrials - Total number of trials
%     orientation - Current orientation angle
%     state - Current state string ('STIMULUS', 'INTERVAL', 'READY')
%     timeElapsed - Time elapsed in current state (seconds)
%     totalProgress - Overall progress (0-1)

% 减少摄像头更新频率
persistent lastCameraUpdate;
if isempty(lastCameraUpdate)
    lastCameraUpdate = 0;
end

currentTime = GetSecs;
if (currentTime - lastCameraUpdate) > 0.1  % 每100ms更新一次摄像头
    % Update camera images
    if isfield(h, 'cam1Available') && h.cam1Available && isvalid(h.cam1)
        try
            img = snapshot(h.cam1);
            set(h.cam1Image, 'CData', img);
        catch
            % If camera errors, do nothing
        end
    end

    if isfield(h, 'cam2Available') && h.cam2Available && isvalid(h.cam2)
        try
            img = snapshot(h.cam2);
            set(h.cam2Image, 'CData', img);
        catch
            % If camera errors, do nothing
        end
    end
    lastCameraUpdate = currentTime;
end

% Update status text - 只显示原始角度
set(h.orientText, 'String', sprintf('Direction: %d°', orientation));
set(h.trialText, 'String', sprintf('Trial: %d / %d', trialNum, totalTrials));
set(h.statusText, 'String', sprintf('Status: %s', state));
set(h.timeText, 'String', sprintf('Time: %d s', timeElapsed));

% Update status color
switch state
    case 'STIMULUS'
        set(h.statusText, 'Color', 'g');
    case 'INTERVAL'
        set(h.statusText, 'Color', 'b');
    otherwise
        set(h.statusText, 'Color', 'k');
end

% Update progress bar
set(h.progressPatch, 'XData', [0 totalProgress totalProgress 0]);

% 使用更高效的刷新方式
drawnow limitrate nocallbacks; % 限制刷新率并跳过回调以提高性能
end