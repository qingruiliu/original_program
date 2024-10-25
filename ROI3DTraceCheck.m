% Loop through the ROI3DWithTraceTable of one session and check the calcium trace of each 3DROI
ROI3DNum = max(ROI3DWithTraceTable.ROI3DIdx);

ROIRegisteredNumber = zeros(ROI3DNum,1);  % the matrix to save the registered status
ROIRegisteredStatus = zeros(ROI3DNum,1);

commonX = 1:5137;

for i = 1:ROI3DNum  % loop through the 3DROI table
    close(gcf)
    tempStruct = ROI3DWithTraceTable.registered_trace_session1(i);
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for j = 1:8  % loop through all the imaging planes
        legendStr = {};
        tempFieldStr = append('Z',num2str(j-1));
        tempValue = tempStruct.(tempFieldStr);
        %tempCount = 0;                                  
        if ~isempty(tempValue)
            %tempCount = tempCount + 1;
            ROIRegisteredNumber(i) = ROIRegisteredNumber(i) + 1;  % if there are registered values
            hold on
            plot(commonX,cell2mat(tempValue)+(80-j*10))
            legendStr{end+1} = tempFieldStr;
        end
        legend(legendStr,'FontSize',16);
    end
    

    titleStr = append('#ROI ',num2str(i));
    title(titleStr,'FontSize',20);
    
    hold off
    
    tempAnswer = questdlg('Is current cell registered?',...
        titleStr,...
        'Not registered', ...
        '1 plane registered', ...
        '>1 planes registered', ...
        'Not registered');
    
    switch tempAnswer
        case 'Not Registered'
            ROIRegisteredStatus(i) = 0;       % 0: no registered
        case '1 plane registered'
            ROIRegisteredStatus(i) = 1;       % 1: only one registered plane
        case '>1 planes registered'
            tempAnswer2 = questdlg('matched or unmatched?', ...
                                    titleStr, ...
                                   'unmatched', ...
                                   'partially matched', ...
                                   'ALL MATCHED', ...
                                   'unmatched');

            switch tempAnswer2
                case 'unmatched'
                    ROIRegisteredStatus(i) = 2;    %2: unmatched planes
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
