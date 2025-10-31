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
    h.cam1 = webcam(1);
    h.cam1Available = true;
catch
    warning('Stimulus monitoring camera (webcam 1) not available');
    h.cam1Available = false;
end

try
    h.cam2 = webcam(2);
    h.cam2Available = true;
catch
    warning('Mouse monitoring camera (webcam 2) not available');
    h.cam2Available = false;
end

% --- 创建监视状态GUI ---
h.monitorFig = figure('Name', 'Experiment Monitor', 'Position', [100, 100, 800, 600], ...
    'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');

% 摄像头预览（如stage2，直接preview，不在主循环中刷新）
h.cam1Panel = uipanel('Parent', h.monitorFig, 'Title', 'Back Camera - Stimulus Monitor', ...
    'Position', [0.025, 0.525, 0.45, 0.45], 'FontSize', 12, 'FontWeight', 'bold');
h.cam2Panel = uipanel('Parent', h.monitorFig, 'Title', 'Front Camera - Mouse Monitor', ...
    'Position', [0.525, 0.525, 0.45, 0.45], 'FontSize', 12, 'FontWeight', 'bold');

h.cam1Axes = axes('Parent', h.cam1Panel, 'Position', [0.05, 0.05, 0.9, 0.9]);
if h.cam1Available
    h.cam1Size = str2double(strsplit(h.cam1.Resolution,'x'));
    h.cam1Image = image(zeros(h.cam1Size),'Parent',h.cam1Axes);
    preview(h.cam1, h.cam1Image); % 只初始化一次，不在主循环中刷新
else
    text(h.cam1Axes, 0.5, 0.5, 'Camera 1 not connected', 'HorizontalAlignment', 'center');
    axis(h.cam1Axes, [0 1 0 1]);
    axis(h.cam1Axes, 'off');
end

h.cam2Axes = axes('Parent', h.cam2Panel, 'Position', [0.05, 0.05, 0.9, 0.9]);
if h.cam2Available
    h.cam2Size = str2double(strsplit(h.cam2.Resolution,'x'));
    h.cam2Image = image(zeros(h.cam2Size),'Parent',h.cam2Axes);
    preview(h.cam2, h.cam2Image); % 只初始化一次，不在主循环中刷新
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

% 移除updateMonitor函数句柄的创建，采用stage2的直接更新方式
% h.updateMonitor = @(trialNum, orientation, state, timeElapsed, totalProgress) updateMonitorGUI(h, trialNum, totalTrials, orientation, state, timeElapsed, totalProgress);

% 预先更新一次GUI - 直接设置，不调用函数
try
    set(h.orientText, 'String', 'Direction: --°');
    set(h.trialText, 'String', 'Trial: 0 / 0');
    set(h.statusText, 'String', 'Status: READY', 'Color', 'k');
    set(h.timeText, 'String', 'Time: 0 s');
    set(h.progressPatch, 'XData', [0 0 0 0]);
catch
    % 忽略GUI更新错误
end
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
    
    % 完全按照stage2的方式：只在试次开始时更新GUI，然后立即执行视觉刺激
    totalProgress = (trialNum - 1) / totalTrials;
    try
        set(h.orientText, 'String', sprintf('Direction: %d°', currentOrientation));
        set(h.trialText, 'String', sprintf('Trial: %d / %d', trialNum, totalTrials));
        set(h.statusText, 'String', 'Status: READY', 'Color', 'k');
        set(h.timeText, 'String', 'Time: 0 s');
        set(h.progressPatch, 'XData', [0 totalProgress totalProgress 0]);
        % 不调用drawnow，让GUI自然刷新，就像stage2一样
    catch
        % 忽略GUI更新错误
    end
    
    % 使用与stage2相同的VBL更新方式
    vbl = Screen('Flip', h.window);
    vblt0 = vbl;
    startTime = vbl;
    
    % 重置相位（与stage2类似）
    h.propertiesMat(1) = 0;
    
    % 刺激呈现循环 - 不做摄像头刷新，最大限度减少延迟
    while vbl - vblt0 <= stimDuration
        Screen('DrawTexture', h.window, h.gratingtex, [], [], displayOrientation, [], [], [], [],...
            kPsychDontDoRotation, h.propertiesMat');
        vbl = Screen('Flip', h.window, vbl + (h.waitframes - 0.5) * h.ifi);
        h.propertiesMat(1) = h.propertiesMat(1) + h.phasePerFrame;
        % 不做drawnow或摄像头刷新
    end
    
    % 刺激结束后简单更新状态，不调用drawnow
    try
        set(h.statusText, 'String', 'Status: STIMULUS Complete', 'Color', 'g');
    catch
        % 忽略GUI更新错误
    end
    
    % Log trial data - 保存相对于实验开始时间的时间戳（秒）
    results.trialLog{trialNum, 1} = trialNum;
    results.trialLog{trialNum, 2} = currentOrientation;
    results.trialLog{trialNum, 3} = round(startTime - experimentStartTime, 4);

    % --- Inter-Stimulus Interval (ISI) ---
    fprintf('ISI period for %f seconds...\n', isiDuration);
    
    % 间隔开始时间
    isiStartTime = Screen('Flip', h.window);
    
    % ISI期间简单更新状态，不调用drawnow
    try
        set(h.statusText, 'String', 'Status: INTERVAL', 'Color', 'b');
    catch
        % 忽略GUI更新错误
    end
    
    % 使用简单的等待，像stage2一样
    while GetSecs < isiStartTime + isiDuration
        WaitSecs(0.1);
    end
    
    % ISI结束后更新进度条，不调用drawnow
    try
        currentProgress = trialNum / totalTrials;
        set(h.progressPatch, 'XData', [0 currentProgress currentProgress 0]);
    catch
        % 忽略GUI更新错误
    end
    
    % 在ISI结束后，利用试次间隙进行数据保存，避免影响时序
    % 每20个试次保存一次数据，减少I/O操作频率
    if mod(trialNum, 20) == 0
        try
            save(results.filename, 'results');
            fprintf('Data saved at trial %d\n', trialNum);
        catch ME
            warning('Failed to save data at trial %d: %s', trialNum, ME.message);
        end
    end
end

% --- End of Experiment ---
try
    set(h.statusText, 'String', 'Status: Experiment Finished', 'Color', 'g');
    set(h.progressPatch, 'XData', [0 1 1 0]);
    drawnow;
catch
    % 忽略GUI更新错误
end

Screen('Flip', h.window);

% 实验结束后关闭摄像头（如stage2）
if h.cam1Available
    stoppreview(h.cam1);
    clear h.cam1;
end
if h.cam2Available
    stoppreview(h.cam2);
    clear h.cam2;
end

% 实验结束后进行最终保存
try
    save(results.filename, 'results');
    fprintf('Final data save completed successfully\n');
catch ME
    warning('Final data save failed: %s', ME.message);
    % 尝试保存到备份文件
    backup_filename = sprintf('backup_%s', results.filename);
    try
        save(backup_filename, 'results');
        fprintf('Data saved to backup file: %s\n', backup_filename);
    catch
        warning('Backup save also failed. Data may be lost.');
    end
end

WaitSecs(2);

% Clean up
sca;
Priority(0);
disp('Experiment finished and data saved.');
fprintf('Results saved to: %s\n', results.filename);