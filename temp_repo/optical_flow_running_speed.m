%% 离线视频分析：选取起始时间 + ROI框选 + 光流法速度分析
clear; clc; close all;

%% 1. 选择视频文件
[fileName, filePath] = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files (*.mp4, *.avi, *.mov)'}, '请选择小鼠行为学视频');
if isequal(fileName, 0)
    disp('取消了文件选择。');
    return;
end
videoFullPath = fullfile(filePath, fileName);
v = VideoReader(videoFullPath);
fps = v.FrameRate;

%% 2. 确定实验正式开始的帧数 (交互式滑块选择)
% 初始化用于存储用户选择结果
startFrameNum = 1; % 默认从第1帧开始

% 创建 UI 图形窗口
figFrame = uifigure('Name', '选择分析起始帧', 'Position', [100, 100, 1000, 650]);
axFrame = uiaxes(figFrame, 'Position', [50, 280, 900, 300]);

% 添加操作提示
uilabel(figFrame, 'Position', [50, 250, 600, 20], ...
    'Text', '控制方式：滑块拖动 或 方向键 (← -1帧, → +1帧, Shift+← -10帧, Shift+→ +10帧)', ...
    'FontColor', [0.4 0.4 0.4], 'FontSize', 11);

% 显示当前帧号和时间
lblFrame = uilabel(figFrame, 'Position', [350, 210, 300, 25], ...
    'Text', 'Current Frame: 1 / ? | Time: 0.00 s', ...
    'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');

% 创建滑块用于选择起始帧
totalFrames = v.NumFrames;
sldFrame = uislider(figFrame, 'Position', [100, 170, 800, 3], ...
    'Limits', [1, max(2, totalFrames)], 'Value', 1);

% 添加帧数输入框 (用户可直接输入帧数)
uilabel(figFrame, 'Position', [50, 130, 100, 20], 'Text', '输入帧数:');
edtFrame = uieditfield(figFrame, 'numeric', 'Position', [150, 125, 80, 28], 'Value', 1, ...
    'Limits', [1, totalFrames]);

% Confirm 按钮
btnConfirm = uibutton(figFrame, 'Text', '确认起始帧', 'FontSize', 13, ...
    'Position', [250, 120, 120, 35], ...
    'BackgroundColor', [0.2, 0.4, 0.8], 'FontColor', 'white');

% Cancel 按钮
btnCancel = uibutton(figFrame, 'Text', '取消', 'FontSize', 13, ...
    'Position', [380, 120, 100, 35], ...
    'BackgroundColor', [0.7, 0.2, 0.2], 'FontColor', 'white');

% 更新显示的回调函数
updateDisplayFcn = @(frameNum) updateFrameDisplay(frameNum, v, axFrame, lblFrame, edtFrame, sldFrame);

% Slider 回调
sldFrame.ValueChangingFcn = @(src, event) updateDisplayFcn(round(event.Value));

% 编辑框回调
edtFrame.ValueChangedFcn = @(src, event) updateDisplayFcn(min(max(round(event.Value), 1), totalFrames));

% 键盘事件处理
figFrame.WindowKeyPressFcn = @(src, event) handleFrameKeyPress(event, sldFrame, totalFrames, updateDisplayFcn);

% Confirm 按钮回调
btnConfirm.ButtonPushedFcn = @(btn, event) set(figFrame, 'UserData', 'confirmed');

% Cancel 按钮回调
btnCancel.ButtonPushedFcn = @(btn, event) set(figFrame, 'UserData', 'cancelled');

% 初始显示
updateDisplayFcn(1);

% 等待用户操作
while isvalid(figFrame)
    userData = get(figFrame, 'UserData');
    if ~isempty(userData) && strcmp(userData, 'confirmed')
        startFrameNum = round(edtFrame.Value);
        delete(figFrame);
        break;
    elseif ~isempty(userData) && strcmp(userData, 'cancelled')
        disp('用户取消了起始帧选择，程序终止。');
        delete(figFrame);
        return;
    end
    pause(0.05);
end

% 根据帧号计算起始时间
startTime = (startFrameNum - 1) / fps;
disp(sprintf('已确认起始帧：第 %d 帧 (时间: %.2f 秒)', startFrameNum, startTime));

%% 3. 交互式选择多个 ROI (感兴趣区域) 及命名 - 手绘模式
% 询问用户要分析多少个ROI
dlg = inputdlg('请输入要分析的ROI个数：', 'ROI个数设置', [1, 30], {'1'});
if isempty(dlg)
    disp('取消了操作，程序终止。');
    return;
end
numROIs = str2double(dlg{1});
if isnan(numROIs) || numROIs < 1
    error('输入的ROI个数无效！');
end

% 初始化ROI存储结构体（存储mask和顶点）
roiData = struct('name', {}, 'mask', {}, 'vertices', {}, 'index', {});

% 将视频跳转到用户指定的起始时间
v.CurrentTime = startTime; 
startFrame = readFrame(v);
startGray = rgb2gray(startFrame);

% 循环让用户手绘每个ROI
for roiIdx = 1:numROIs
    % 创建窗口用于手绘
    hFig = figure('Name', sprintf('手绘ROI %d/%d', roiIdx, numROIs), 'NumberTitle', 'off', 'Position', [100, 100, 900, 700]);
    imshow(startGray);
    title(sprintf('ROI %d: 用鼠标绘制感兴趣区域 (绘制完成后双击或按Enter确认)', roiIdx), 'FontSize', 12);
    
    % 使用 drawfreehand 让用户手绘自由形状
    roi = drawfreehand('Color', 'r', 'LineWidth', 2);
    
    % 获取绘制的多边形顶点
    vertices = roi.Position;
    
    % 创建mask
    [height, width] = size(startGray);
    mask = poly2mask(vertices(:, 1), vertices(:, 2), height, width);
    
    close(hFig);
    
    % 弹出对话框询问ROI名称
    roiName = inputdlg(sprintf('请输入第 %d 个ROI的名称 (例如：左轮、右轮、舔水器等)：', roiIdx), ...
                       '命名ROI', [1, 50], {sprintf('ROI_%d', roiIdx)});
    if isempty(roiName)
        roiName = {sprintf('ROI_%d', roiIdx)}; % 如果用户取消，使用默认名称
    end
    
    % 存储ROI信息
    roiData(roiIdx).name = roiName{1};
    roiData(roiIdx).mask = mask;
    roiData(roiIdx).vertices = vertices;
    roiData(roiIdx).index = roiIdx;
    
    disp(['ROI ', num2str(roiIdx), ' - 名称: "', roiData(roiIdx).name, '" - 已绘制']);
end

%% 4. 初始化光流对象与变量
% 为每个ROI初始化光流对象（使用元胞数组，因为 opticalFlowFarneback 对象不支持普通数组索引赋值）
opticFlows = cell(numROIs, 1);
for i = 1:numROIs
    opticFlows{i} = opticalFlowFarneback;
end

% 【极其关键的一步】：将视频时间再次重置回用户选择的起始点
% 这样接下来的 while 循环才是真正从你选定的那一帧开始分析
v.CurrentTime = startTime; 

% 预估剩余帧数以分配内存（提高 MATLAB 运行速度）
remainingDuration = v.Duration - startTime;
estFrames = ceil(remainingDuration * fps);
% 为每个ROI分别存储运动强度数据
motionIntensities = zeros(estFrames, numROIs);

%% 5. 循环处理每一帧
disp('开始计算光流，请耐心等待...');
hWait = waitbar(0, '正在分析视频，请稍候...');
frameCount = 0;

while hasFrame(v)
    frameCount = frameCount + 1;
    frame = readFrame(v);
    gray = rgb2gray(frame);
    
    % 对每个ROI分别计算光流（使用mask而非imcrop）
    for roiIdx = 1:numROIs
        flow = estimateFlow(opticFlows{roiIdx}, gray);
        
        % 应用mask：只在ROI区域内计算运动强度
        maskedMagnitude = flow.Magnitude;
        maskedMagnitude(~roiData(roiIdx).mask) = 0;  % 掩码外设为0
        
        % 计算mask内的平均运动强度
        validPixels = maskedMagnitude(roiData(roiIdx).mask);
        if ~isempty(validPixels)
            motionIntensities(frameCount, roiIdx) = mean(validPixels);
        else
            motionIntensities(frameCount, roiIdx) = 0;
        end
    end
    
    if mod(frameCount, 100) == 0
        % 更新进度条（使用预估总帧数）
        progress = min(frameCount / estFrames, 1); 
        waitbar(progress, hWait, sprintf('处理进度: %d 帧', frameCount));
    end
end
close(hWait);

% 截断多余的预分配内存（防止预估帧数与实际提取帧数有细微偏差）
motionIntensities = motionIntensities(1:frameCount, :);

%% 6. 数据平滑、时间轴对齐与绘图
% 创建以 0 为起点的时间轴，方便与你的 behavioral program 的 tic/toc 对齐
timeAxis = (0 : frameCount-1) / fps; 

% 0.5秒滑动平均平滑去噪
smoothWindow = round(fps / 2);
smoothSpeeds = zeros(frameCount, numROIs);
for roiIdx = 1:numROIs
    smoothSpeeds(:, roiIdx) = movmean(motionIntensities(:, roiIdx), smoothWindow);
end

% 定义颜色用于不同的ROI
colors = lines(numROIs); % 自动生成numROIs种颜色

% 绘图展示
fig = figure('Name', 'Multi-ROI Running Speed Analysis', 'Color', 'w');
hold on;
legendLabels = {};
for roiIdx = 1:numROIs
    plot(timeAxis, smoothSpeeds(:, roiIdx), 'LineWidth', 1.5, 'Color', colors(roiIdx, :));
    legendLabels{roiIdx} = roiData(roiIdx).name;
end
hold off;
legend(legendLabels, 'Location', 'best', 'FontSize', 10);
xlabel('Time from experiment start (seconds)', 'FontSize', 12);
ylabel('Motion Intensity (pixels/frame)', 'FontSize', 12);
title(['Mouse Running Activity - ', fileName], 'FontSize', 14, 'Interpreter', 'none');
grid on; box off;

%% 7. 导出数据
[~, name, ~] = fileparts(fileName);
saveNameMat = fullfile(filePath, [name, '_MultiROI_SpeedData.mat']);
saveNameCsv = fullfile(filePath, [name, '_MultiROI_SpeedData.csv']);

% 保存 .mat 格式供后续与 mLatency 联合分析
% 同时保存ROI信息方便后续追踪
save(saveNameMat, 'timeAxis', 'motionIntensities', 'smoothSpeeds', 'roiData', 'fps', 'startTime');

% 构建表格并导出 .csv 格式
% 创建一个包含所有ROI数据的表格
tableData = table(timeAxis');
tableData.Properties.VariableNames = {'Time_sec'};

for roiIdx = 1:numROIs
    % 原始运动数据列
    colNameRaw = ['Raw_Motion_' strrep(roiData(roiIdx).name, ' ', '_')];
    % 平滑运动数据列
    colNameSmooth = ['Smoothed_Motion_' strrep(roiData(roiIdx).name, ' ', '_')];
    
    tableData.(colNameRaw) = motionIntensities(:, roiIdx);
    tableData.(colNameSmooth) = smoothSpeeds(:, roiIdx);
end

writetable(tableData, saveNameCsv);

% 也保存一份ROI配置文件便于查阅
saveNameConfig = fullfile(filePath, [name, '_ROI_Config.txt']);
fid = fopen(saveNameConfig, 'w');
fprintf(fid, '========== ROI 配置信息 (手绘多边形) ==========\n');
fprintf(fid, '视频文件: %s\n', fileName);
fprintf(fid, '起始时间: %.2f 秒\n', startTime);
fprintf(fid, '帧率: %.2f fps\n\n', fps);
fprintf(fid, 'ROI 详细信息:\n');
for roiIdx = 1:numROIs
    fprintf(fid, '\nROI %d: %s\n', roiIdx, roiData(roiIdx).name);
    fprintf(fid, '多边形顶点坐标 (x, y):\n');
    vertices = roiData(roiIdx).vertices;
    for vIdx = 1:size(vertices, 1)
        fprintf(fid, '  顶点 %d: (%.1f, %.1f)\n', vIdx, vertices(vIdx, 1), vertices(vIdx, 2));
    end
end
fclose(fid);

disp(['分析完成！起始点已校准。']);
disp(['数据已保存至:']); 
disp(['  MAT 文件: ', saveNameMat]);
disp(['  CSV 文件: ', saveNameCsv]);
disp(['  配置文件: ', saveNameConfig]);

%% ========== 辅助函数 ==========

% 更新帧显示的回调函数
function updateFrameDisplay(frameNum, vObj, axObj, lblObj, edtObj, sldObj)
    frameNum = max(1, min(round(frameNum), vObj.NumFrames));
    img = read(vObj, frameNum);
    imshow(img, 'Parent', axObj);
    
    % 计算时间
    frameRate = vObj.FrameRate;
    timeInSec = (frameNum - 1) / frameRate;
    
    % 更新标签
    lblObj.Text = sprintf('Current Frame: %d / %d | Time: %.2f s', frameNum, vObj.NumFrames, timeInSec);
    
    % 同步更新滑块和编辑框（避免循环触发）
    if abs(sldObj.Value - frameNum) > 0.1
        sldObj.Value = frameNum;
    end
    if edtObj.Value ~= frameNum
        edtObj.Value = frameNum;
    end
end

% 处理键盘事件
function handleFrameKeyPress(event, sldObj, maxFrames, refreshFcn)
    step = 1;
    if any(strcmp(event.Modifier, 'shift'))
        step = 10;
    end
    
    switch event.Key
        case 'leftarrow'
            newValue = max(1, sldObj.Value - step);
            sldObj.Value = newValue;
            refreshFcn(newValue);
        case 'rightarrow'
            newValue = min(maxFrames, sldObj.Value + step);
            sldObj.Value = newValue;
            refreshFcn(newValue);
    end
end