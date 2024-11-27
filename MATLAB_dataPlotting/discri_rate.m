%% correct rate during discrimination training
correctRateMat = NaN(11,14);

M3 = [146 153 154 241 243 289];
M4 = [154 154 144 167 167 178 172 176 212 216 256 247 244];
M7 = [173 199 233 273 285 293];
M8 = [155 164 218 254 250 270];
M9 = [139 147 152 157 159 175 207 194 248 261 290];
M12 = [169 185 218 285 289 290];
M15 = [160 252 283 281];
M17 = [155 154 173 208 252 264 280];
M18 = [168 171 183 222 212 236 233 234 248 273 261];
M20 = [140 161 155 155 152 175 190 193 190 180 253 242 274];
M21 = [154 156 149 163 173 201 263 291 267]; 

correctRateMat(1,1:length(M3)) = M3;
correctRateMat(2,1:length(M4)) = M4;
correctRateMat(3,1:length(M7)) = M7;
correctRateMat(4,1:length(M8)) = M8;
correctRateMat(5,1:length(M9)) = M9;
correctRateMat(6,1:length(M12)) = M12;
correctRateMat(7,1:length(M15)) = M15;
correctRateMat(8,1:length(M17)) = M17;
correctRateMat(9,1:length(M18)) = M18;
correctRateMat(10,1:length(M20)) = M20;
correctRateMat(11,1:length(M21)) = M21;

correctRateMat = correctRateMat/300 * 100; % convert to percentage

% plot correct rate
figure;
rateAxes = axes;
commonX = 1:14;
for i = 1:11
    plot(rateAxes,commonX,correctRateMat(i,:),'-','lineWidth',2,'Color',[0.5 0.5 0.5 0.4]);
    hold on
end
plot(rateAxes,commonX,correctRateMat(6,:),'-','lineWidth',2,'Color',[0 0 0]); %example data
hold off
%set the tick label length as 0
set(rateAxes,'TickLength',[0 0])
xlabel(rateAxes,'Training day');
xticks([1 5 10 14]);
xticklabels({'1','5','10','14'});
xlim([1 14])

ylabel(rateAxes,'Correct rate (%)');
yticks([50 80 100])
ylim([45 100])
yline(80,'--','LineWidth',1); %expert line
set(rateAxes,'FontSize',18)