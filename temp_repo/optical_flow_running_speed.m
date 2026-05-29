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

%% 2. 确定实验正式开始的时间点
% 调出 MATLAB 自带的视频播放器，供用户寻找起始时间
disp('正在打开视频播放器，请寻找实验开始的准确时间（秒）...');
implayHandle = implay(videoFullPath);

% 弹出输入对话框
prompt = {'请在播放器中查看，并输入实验正式开始的时间（单位：秒）：', '找到时间后输入，您可以直接关闭视频播放器。'};
dlgtitle = '设置分析起始点';
dims = [1 50; 1 50];
definput = {'0', ''};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('取消了输入，程序终止。');
    if isvalid(implayHandle), close(implayHandle); end
    return;
end

startTime = str2double(answer{1});
if isnan(startTime) || startTime < 0 || startTime > v.Duration
    error('输入的时间无效，请重新运行脚本。');
end

% 如果播放器还没关，程序自动帮你关掉以释放内存
if isvalid(implayHandle)
    close(implayHandle);
end

%% 3. 交互式选择多个 ROI (感兴趣区域) 及命名
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

% 初始化ROI存储结构体
roiData = struct('name', {}, 'rect', {}, 'index', {});

% 将视频跳转到用户指定的起始时间
v.CurrentTime = startTime; 
startFrame = readFrame(v);
startGray = rgb2gray(startFrame);

% 循环框选每个ROI
for roiIdx = 1:numROIs
    % 弹出窗口让用户框选区域
    hFig = figure('Name', sprintf('起始帧 (%.1f 秒) - ROI %d/%d', startTime, roiIdx, numROIs), 'NumberTitle', 'off');
    imshow(startGray);
    title(sprintf('请框选第 %d 个ROI，然后【双击框内】确认', roiIdx));
    [~, roiRect] = imcrop(hFig); 
    close(hFig);
    
    roiRect = round(roiRect);
    
    % 弹出对话框询问ROI名称
    roiName = inputdlg(sprintf('请输入第 %d 个ROI的名称 (例如：左轮、右轮、中心等)：', roiIdx), ...
                       '命名ROI', [1, 50], {sprintf('ROI_%d', roiIdx)});
    if isempty(roiName)
        roiName = {sprintf('ROI_%d', roiIdx)}; % 如果用户取消，使用默认名称
    end
    
    % 存储ROI信息
    roiData(roiIdx).name = roiName{1};
    roiData(roiIdx).rect = roiRect;
    roiData(roiIdx).index = roiIdx;
    
    disp(['ROI ', num2str(roiIdx), ' - 名称: "', roiData(roiIdx).name, ...n'" - 坐标: X=', num2str(roiRect(1)), ', Y=', num2str(roiRect(2)), ...
          ', 宽=', num2str(roiRect(3)), ', 高=', num2str(roiRect(4))]);
end

%% 4. 初始化光流对象与变量
% 为每个ROI初始化光流对象
opticFlows = repmat(opticalFlowFarneback, numROIs, 1);

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
    
    % 对每个ROI分别裁剪并计算光流
    for roiIdx = 1:numROIs
        croppedGray = imcrop(gray, roiData(roiIdx).rect);
        flow = estimateFlow(opticFlows(roiIdx), croppedGray);
        % 记录当前ROI的运动强度
        motionIntensities(frameCount, roiIdx) = mean(flow.Magnitude(:));
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
fprintf(fid, '========== ROI 配置信息 ==========\n');
fprintf(fid, '视频文件: %s\n', fileName);
fprintf(fid, '起始时间: %.2f 秒\n', startTime);
fprintf(fid, '帧率: %.2f fps\n\n', fps);
fprintf(fid, 'ROI 详细信息:\n');
for roiIdx = 1:numROIs
    fprintf(fid, '\nROI %d: %s\n', roiIdx, roiData(roiIdx).name);
    fprintf(fid, '  坐标 X: %d\n', roiData(roiIdx).rect(1));
    fprintf(fid, '  坐标 Y: %d\n', roiData(roiIdx).rect(2));
    fprintf(fid, '  宽度: %d\n', roiData(roiIdx).rect(3));
    fprintf(fid, '  高度: %d\n', roiData(roiIdx).rect(4));
end
fclose(fid);

disp(['分析完成！起始点已校准。']);
disp(['数据已保存至:']); 
disp(['  MAT 文件: ', saveNameMat]);
disp(['  CSV 文件: ', saveNameCsv]);
disp(['  配置文件: ', saveNameConfig]);