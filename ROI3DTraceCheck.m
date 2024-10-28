% Loop through the ROI3DWithTraceTable of one session and check the calcium trace of each 3DROI
ROI3DNum = max(ROI3DWithTraceTable.ROI3DIdx);

ROIRegisteredNumber = zeros(ROI3DNum,1);  % the matrix to save the registered status
ROIRegisteredStatus = zeros(ROI3DNum,1);

commonX = 1:1000;    %only plot the first 1000 frames

for i = 1:ROI3DNum  % loop through the 3DROI table
    close(gcf)
    tempStruct = ROI3DWithTraceTable.registered_trace_session1(i);
    figure('units','normalized','outerposition',[0 0.5 0.5 0.5])    %show the new figure at left up corner
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
    ylabel('Deconv F (a.u.)');
    timeScale = 20 * 2.2;     %scale bar of time: 10 seconds * frequency
    ampScale  = 5;     %scale bar of amplitude: 5 a.u.

    xPos = commonX(end) - timeScale - 5;  %start point of time scale
    yPos = ampScale;                      %amp scale start from 0

    plot([xPos, xPos + timeScale],[yPos, yPos],'k','LineWidth',2);
    text(xPos + timeScale/2, yPos - 2, [num2str(timeScale/2.2),' s'],...
         'HorizontalAlignment','center');       %plot time scale bar

    plot([xPos + timeScale, xPos + timeScale],[yPos, yPos + ampScale],'k','LineWidth',2);
    text(xPos + timeScale + 2, yPos + ampScale/2,[num2str(ampScale),' a.u.'],...
         'HorizontalAlignment','left','Rotation',90);
    axis off

    legend(legendStr,'FontSize',16);
    titleStr = append('#ROI ',num2str(i));
    title(titleStr,'FontSize',20);
    
    hold off
    
    if tempCount == 0
        ROIRegisteredStatus(i) = 0;    %not registered to any planes

    elseif tempCount == 1   
        ROIRegisteredStatus(i) = 1;    %registered to one plane
        
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
