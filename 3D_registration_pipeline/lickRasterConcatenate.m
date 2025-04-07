%This program is used for the output of the behavioral data. 
% 2024.3.28 Liu @UTokyo

%% extract the trial raster in the .mat file

% Load the selected .mat file and get the name of it
fileName = uigetfile('*.mat','MultiSelect','on');
wb = waitbar(0,'Loading data...');
for i = 1:length(fileName)
    waitbar(i/length(fileName),wb,'Loading data...');
    tempStr = fileName{i};
    load(tempStr)
    tempFileName = tempStr(1:end-4); 

    axesPosition = get(h.trialRaster, 'Position');

    % Create a new figure and resize, remember to open the UI graph
    hNewFigure = figure;
    newRasterName = append(tempFileName,'Raster');
    set(hNewFigure,'Name',newRasterName);

    hNewAxes = copyobj(h.trialRaster,hNewFigure);

    % Calculate the position to center the axes in the new figure
    newAxesPosition = [0.1, 0.1, 0.8, 0.8];

    % Copy the axes to the new figure
    set(hNewAxes, 'Position', newAxesPosition);
    titleStr = [newRasterName, '.fig'];
    saveas(hNewAxes, titleStr)
end
close(wb)

%% Prompt user to select multiple files for batch processing
[fileNames, pathName] = uigetfile('*.fig', 'Select raster files', 'MultiSelect', 'on');

if ischar(fileNames)
    fileNames = {fileNames}; % Ensure fileNames is a cell array
end

% Y-offset for each day's raster plot
yOffsets = 300 * (0:numel(fileNames)-1);

fig = figure;
hNewAxes = axes(fig,"Position",[0.15 0.25 0.48 0.6]);
% Open the first figure and copy its plot information
for j = 1:length(fileNames)
    hFig1 = openfig(fileNames{j},'invisible');
    originalAxes = findobj(hFig1,'Type','axes');
    originalLines = findobj(originalAxes,'Type','line');


  for i = 1:numel(originalLines)
   xdata = get(originalLines(i),'XData');
   ydata = get(originalLines(i),'YData')+ yOffsets(j);

   plot(hNewAxes, xdata, ydata,...
       'Color',get(originalLines(i),'Color'),...
       'LineStyle',get(originalLines(i),'LineStyle'), ...
       'LineWidth',get(originalLines(i),'LineWidth'), ...
       'Marker', get(originalLines(i),'Marker'),...
       'MarkerSize', get(originalLines(i),'MarkerSize'));
   hold(hNewAxes,'on');
  end
end

%% Loop through each file and copy the plot information
%
set(gca,'YDir','reverse')
xlim([-1 9])
xticks([0 1 5 9])
ylim([0 1800])
yticks(0:300:1800)
set(gca,'TickLength',[0.005, 0.005])
xlabel('Time (s)')
ylabel('Trial number')
%title('trial raster of C3M2')
xline(1,'-','Color',[0 0 0],'LineWidth',1)
xline(5,'-','Color',[0 0 0],'LineWidth',1)
ha = xregion([0 1],'FaceColor',[0.5 0.5 0.5]);
alpha(ha,0.2)
set(gca,'TickLength',[0 0])
set(gca,'FontSize',16)
for i = 1:6
    yline(i*300,'--','LineWidth',0.5)
end
%% add correct rate plot on the right (across sessions)
rateAxes = axes(fig,'Position',[0.67 0.25 0.2 0.6],...
                'XDir','reverse');
correctRate = [169 185 218 285 289 290];
correctRate = correctRate/300 * 100;
correctRateY =0.5:1:5.5;
plot(rateAxes,correctRate,correctRateY,'k.-','lineWidth',1.5,'MarkerSize',15);
set(rateAxes,'YDir','reverse');
set(rateAxes,'YAxisLocation','Right')
xlabel(rateAxes,'Correct rate (%)');
ylabel(rateAxes,'Training day');
xline(80,'--','LineWidth',0.5);
ylim([0 6]);
yticks(0.5:1:5.5);
yticklabels({'1','2','3','4','5','6'})
xlim([50 100]);
xticks([50 80 100]);
set(gca,'TickLength',[0 0])
set(gca,'FontSize',16)

set(gcf,'color','none')

