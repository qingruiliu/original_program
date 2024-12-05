%% select the result data at same contrast leve;
%hit
disp('----------choose the Hit----------')
[fileName,path] = uigetfile('*Hit*.mat');
cd(path)
load(fileName)
Hit_temp = tempVar;

%Miss
disp('----------choose the Miss----------')
[fileName,path] = uigetfile('*Miss*.mat');
cd(path)
load(fileName)
Miss_temp = tempVar;

%FA
disp('----------choose the FA----------')
[fileName,path] = uigetfile('*FA*.mat');
cd(path)
load(fileName)
FA_temp = tempVar;

%CR
disp('----------choose the CR----------')
[fileName,path] = uigetfile('*CR*.mat');
cd(path)
load(fileName)
CR_temp = tempVar;

%%