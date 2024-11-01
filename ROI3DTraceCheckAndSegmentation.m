%% Loop through the ROI3DWithTraceTable of one session and check the calcium trace of each 3DROI
ROI3DNum = max(ROI3DWithTraceTable.ROI3DIdx);

ROIRegisteredNumber = zeros(ROI3DNum,1);  % the matrix to save the registered status
ROIRegisteredStatus = zeros(ROI3DNum,1);

sessionNum = questdlg('Which session of data you want to check:','session check',...
                                   'session1','session2','session3','session1');      %choose the session to check
  
sessionStr = append('registered_trace_',sessionNum);                                  %create the session string

timeLagM3 = [2.666666 0.9 4.1666666];      %the imaging-behavioral lapse of 3 sessions of M3

tempTimeLag = timeLagM3(str2num(extractAfter(sessionNum,'n')));   %after selection, the lag of current session is selected.
%% load the timeStamp and behavioral output for segmentation
disp('--------select ALL the timestamp .txt files of the current session-------------');
[timeStamps, timeStampPath] = uigetfile('*.txt',  'All Files (*.*)','MultiSelect','on');
cd(timeStampPath);

timeStampTable = table();                       %create a new timeStampTable of current session
for i = 1:length(filename)
    tempStampStr = append('Z',num2str(i-1));
    tempData = importdata(filename{i});         %the first row of time stamp information is necessary
    tempData(:,2) = [];
    
    %make each timestamp evenly spaced
    frameTime = (max(tempData)-tempData(1))/(length(tempData)-1);
    tempData = (tempData(1):frameTime:max(tempData))';
    timeStampTable.(tempStampStr) = tempData;
end

%% load the behavioral data 
disp('---------select the behavioral data of current session-----------------');
[bhfile,bhpath] = uigetfile('.mat','Select the behavioral data of current session');
cd(bhpath);
load(befile);

%% prepare the plot for each ROI's registered traces

commonX = 1:1000;    %only plot the first 1000 frames

for i = 91:ROI3DNum  % loop through the 3DROI table
    close(gcf)
    tempStruct = ROI3DWithTraceTable.(sessionStr)(i);       %check registered session2 trace 24.10.31
    f = figure('units','normalized','outerposition',[0 0.5 0.5 0.5]);    %show the new figure at left up corner
    legendStr = {};                                                 %pre assign the empty legend string
    tempCount = 0;          %counter of the imaging planes with calcium trace

    for j = 1:8  % loop through all the imaging planes
        tempFieldStr = append('Z',num2str(j-1));
        tempValue = cell2mat(tempStruct.(tempFieldStr));                                
        if ~isempty(tempValue)
            tempCount = tempCount + 1;
            ROIRegisteredNumber(i) = ROIRegisteredNumber(i) + 1;  % if there are registered values
            hold on
            plot(commonX,tempValue(1:1000)+(80-j*10))
            legendStr{end+1} = tempFieldStr;                      % add the legend string
        end
    end
    
    %add the scale bar and change the info of axes
    xlabel('Time (s)');
    xlim([0 1000]);
    ylabel('Deconv F (a.u.)');
    ylim([-10 inf])

    timeScale = 20 * 2.2;     %scale bar of time: 10 seconds * frequency
    ampScale  = 10;     %scale bar of amplitude: 5 a.u.

    xPos = commonX(end) - timeScale - 5;  %start point of time scale
    yPos = ampScale;                      %amp scale start from 0

    plot([xPos, xPos + timeScale],[yPos-10, yPos-10],'k','LineWidth',1);
    text(xPos + timeScale/2, yPos - 12, [num2str(timeScale/2.2),' s'],...
         'HorizontalAlignment','center');       %plot time scale bar

    plot([xPos + timeScale, xPos + timeScale],[yPos - 10, yPos - 10 + ampScale],'k','LineWidth',1);
    text(xPos + timeScale + 10, yPos - 10 + ampScale/2,[num2str(ampScale),' a.u.'],...
         'HorizontalAlignment','center','Rotation',90);
    axis off


    legend(legendStr,'FontSize',12,'Location','best');       %avoid blocking the plotting area
    titleStr = append('ROI #',num2str(i));
    title(titleStr,'FontSize',20);
    
    hold off

    saveStr = append(titleStr,'.png');
    saveas(f,saveStr)
    if tempCount == 0
        ROIRegisteredStatus(i) = 0;    %not registered to any planes
        ROI3DWithTraceTable.(sessionStr)(i).selected = [];   %11.1 add the selected field to the struct,when segmentation, only use the .selected field
    elseif tempCount == 1   
        ROIRegisteredStatus(i) = 1;    %registered to one plane
        tempTrace = ROI3DWithTraceTable.(sessionStr)(i).(legendStr{1});         %give the .selected field with the only field with value

        ROI3DWithTraceTable.(sessionStr)(i).selected = selectedTraceSegment(tempTrace,legendStr{1});  %use the function to segment the selected trace
    elseif tempCount >= 2
        tempAnswer2 = questdlg('matched or unmatched?', ...
                                titleStr, ...
                                'unmatched', ...
                                'partially matched', ...
                                'ALL MATCHED', ...
                                'unmatched');
        switch tempAnswer2
            case 'unmatched'
                ROIRegisteredStatus(i) = 2;    %2: unmatched multiple planes
            case 'partially matched'
                ROIRegisteredStatus(i) = 3;    %3: partially matched planes
            case 'ALL MATCHED'
                ROIRegisteredStatus(i) = 4;    %4: all matched planes
        end
        
        [indx,tf] = listdlg('PromptString',{'Select the plane represent this 3D-ROI'},...
                            'ListString',legendStr,'SelectionMode','single',...
                            'ListSize',[500 300]);     %use the fields with value to create a list dialog box
        
        if tf 
            tempTrace = ROI3DWithTraceTable.(sessionStr)(i).(legendStr{indx});
            ROI3DWithTraceTable.(sessionStr)(i).selected = selectedTraceSegment(tempTrace,legendStr{indx}); %use the function to segment the selected trace
        end
        
    end

    tempAnswer3 = questdlg('Move to the next ROI?', ...
                            'NEXT', ...
                            'Next ROI', ...
                            'Stop the program', ...
                            'Next ROI');
    switch tempAnswer3
        case 'Next ROI'
            continue
        case 'Stop the program'
            break
    end
    
end


function segmentedTrace = selectedTraceSegment(trace,depth)


end