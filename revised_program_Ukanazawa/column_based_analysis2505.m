%% selected the trial-segmented table including transformed ROI coordinates 
waitfor(msgbox('Select the trial-segmented table including transformed ROI coordinates'))
[filename, pathname] = uigetfile({'*.mat'},'Select the trial-segmented table including transformed ROI coordinates');
cd(pathname);
load(filename);

trans = ROISegTraceTable.TransformedCenter;
transCenter = trans(:,1:2);
%% autocorrelograms analysis

% Load your neuron positions: Nx2 matrix
% Assume transCenter = [x, y];
% Load or define this before running the code.

% Parameters
bin_width = 10;           % Width of each bin in micrometers
max_distance = 300;       % Max distance to evaluate in micrometers
orientation_deg = 0;      % Orientation of periodicity in degrees (0 = horizontal)

% Convert orientation to radians
theta = deg2rad(orientation_deg);

% Orientation vectors
dir_vector = [cos(theta), sin(theta)];       % Along periodicity
perp_vector = [-sin(theta), cos(theta)];     % Perpendicular to periodicity

% Bin centers along direction vector
distances = -max_distance:bin_width:max_distance;
cell_density = zeros(size(distances));

% Iterate over all neurons as reference points
num_cells = size(transCenter, 1);
for i = 1:num_cells
    ref_cell = transCenter(i, :);

    % Exclude the reference cell
    other_cells = transCenter;
    other_cells(i, :) = [];  % Remove reference cell
    rel_positions = other_cells - ref_cell;


    % Project relative positions onto direction and perpendicular axes
    parallel_proj = rel_positions * dir_vector';
    perp_proj = abs(rel_positions * perp_vector');  % Distance to the bin line

    % Loop over distance bins
    for j = 1:length(distances)
        d = distances(j);
        
        % Define bin region: band at position d along dir_vector
        bin_center = d;
        in_bin = (abs(parallel_proj - bin_center) < bin_width/2) & (perp_proj < bin_width/2);
        
        % Count how many other cells fall in the bin
        cell_density(j) = cell_density(j) + sum(in_bin);
    end
end

% Normalize by number of reference cells and bin width (gives density)
cell_density = cell_density / (num_cells * bin_width);  % units: cells/μm

% Plot as bar/step (each bin为一条水平线)
figure;
for j = 1:length(distances)-1
    x = [distances(j), distances(j+1)];
    y = [cell_density(j), cell_density(j)];
    plot(x, y, 'k-', 'LineWidth', 2); hold on;
end
xlabel('Distance along periodicity (\mum)');
ylabel('Cell density (cells/\mum)');
title('Autocorrelogram');
grid on;
hold off;

% 计算XY平面距离并聚类
xyCoord = transCenter(:,1:2);
numCells = size(xyCoord,1);

% 计算距离矩阵
D = squareform(pdist(xyCoord));

% 在矩阵D的每一列中寻找距离小于20的元素的索引
closeIdxCell = cell(1, numCells);
for i = 1:numCells
    closeIdxCell{i} = find(D(:,i) < 20 & D(:,i) > 0); % 排除自身
end
% closeIdxCell{i} 存储第i个细胞距离小于20的所有细胞索引

% 距离小于20的视为同一微柱，使用连通分量聚类
G = D < 20 & D > 0; % 排除自身
G = sparse(G);
[~, microcolumnIdx] = graphconncomp(G, 'Weak', true);

% microcolumnIdx为每个细胞的微柱编号
% 可视化时用不同颜色区分微柱
figure;
scatter3(transCenter(:,1),transCenter(:,2),transCenter(:,3),30,microcolumnIdx,'filled','MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
colormap(jet(max(microcolumnIdx)));
colorbar;
title('Microcolumn Clustering by XY Distance');