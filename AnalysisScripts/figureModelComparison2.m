close all;
axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
AMmarker = 'd';
IOmarker = '*';
markerEdgeColors = [0 0 0.35; 1 0.4 0];
markerLineWidth = 1;
markerSize = 6;

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);
propL = ADT.L./(ADT.L + ADT.M);
minPropL = ADT.L./(ADT.L + ADT.M + ADT.X);
maxPropL = (ADT.L + ADT.X)./(ADT.L + ADT.M + ADT.X);

S = [380 1 401];
wvl = SToWls(S);

% get real data
convertThresholdsToQuanta;
realData = retIrradGreenQuantaPerSecM2./retIrradRedQuantaPerSecM2;

% Import ideal observer data
% IO555 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO555.mat');
% IO559 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO559.mat');
% IO563 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO563.mat');
% 
% IO_RedThresh555 = 10.^IO555.logThreshold680;
% IO_GreenThresh555 = 10.^IO555.logThreshold543;
% IO_RedThresh559 = 10.^IO559.logThreshold680;
% IO_GreenThresh559 = 10.^IO559.logThreshold543;
% IO_RedThresh563 = 10.^IO563.logThreshold680;
% IO_GreenThresh563 = 10.^IO563.logThreshold543;
% 
% IOSensitivityRatio555 = IO_GreenThresh555./IO_RedThresh555;
% IOSensitivityRatio559 = IO_GreenThresh559./IO_RedThresh559;
% IOSensitivityRatio563 = IO_GreenThresh563./IO_RedThresh563;

%additive model

[RedSensitivity555, GreenSensitivity555, ~] = simplePigmentModel(555.5);
[RedSensitivity559, GreenSensitivity559, ~] = simplePigmentModel;
[RedSensitivity563, GreenSensitivity563, ~] = simplePigmentModel(563.4);

propLFine = linspace(0, 1, 1e3);


additiveSensitivityRatio555 = interp1(propLFine, RedSensitivity555, propL)./interp1(propLFine, GreenSensitivity555, propL);
additiveSensitivityRatio559 = interp1(propLFine, RedSensitivity559, propL)./interp1(propLFine, GreenSensitivity559, propL);
additiveSensitivityRatio563 = interp1(propLFine, RedSensitivity563, propL)./interp1(propLFine, GreenSensitivity563, propL);

maxSensitivityRatio = max([realData; IOSensitivityRatio555; IOSensitivityRatio559; IOSensitivityRatio563; additiveSensitivityRatio555; additiveSensitivityRatio559; additiveSensitivityRatio563]);
minSensitivityRatio = min([realData; IOSensitivityRatio555; IOSensitivityRatio559; IOSensitivityRatio563; additiveSensitivityRatio555; additiveSensitivityRatio559; additiveSensitivityRatio563]);

L_lambdaMax = [563.4 558.9];

figure;
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5    13.5   13.5], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters');


ax1 = axes('Position', [0.125 0.4, 0.35 0.35], 'LineWidth', defaultLineWidth, 'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
xlabel(ax1, 'Obs. sensitivity ratio', 'FontSize', labelFontSize);
ylabel(ax1, 'Pred. sensitivity ratio', 'FontSize', labelFontSize);
xlim(ax1, [0 0.035]); ylim(ax1, [0 0.035]);


% ax2 = axes('Position', [0.625, 0.6, 0.35 0.35], 'LineWidth', defaultLineWidth, 'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
% xlabel(ax2, 'Proportion L', 'FontSize', labelFontSize);
% ylabel(ax2, 'Residual', 'FontSize', labelFontSize);
% xlim(ax2, [0 1.05]); ylim(ax2, [-0.015 0.015]);
% 
% 
% ax3 = axes('Position', [0.125, 0.125, 0.35 0.35], 'LineWidth', defaultLineWidth,  'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
% xlabel(ax3, 'A.M. sensitivity ratio', 'FontSize', labelFontSize);
% ylabel(ax3, 'I.O. sensitivity ratio', 'FontSize', labelFontSize);
% xlim(ax3, [0 0.035]); ylim(ax3, [0 0.035]);

ax4 = axes('Position', [0.625, 0.4, 0.35 0.35], 'LineWidth', defaultLineWidth,  'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
xlabel(ax4, '\lambda_{max} (nm)', 'FontSize', labelFontSize);
ylabel(ax4, 'RMSE', 'FontSize', labelFontSize);
xlim(ax4, [550 570]);

hold(ax1, 'on');
% hold(ax2, 'on');
% hold(ax3, 'on');
hold(ax4, 'on');






for subj = 1:numel(subjects)

% create figure for each subject
subjIdx = strcmpi(ADT.SubjectID, subjects{subj});


%% Set up axes



switch subjects{subj}
    case '10001R'
        IOSensitivityRatioBestLambdaMax = IOSensitivityRatio563(subjIdx);
        additiveSensitivityRatioBestLambdaMax = additiveSensitivityRatio563(subjIdx);
        
    case '20217R'
        IOSensitivityRatioBestLambdaMax = IOSensitivityRatio559(subjIdx);
        additiveSensitivityRatioBestLambdaMax = additiveSensitivityRatio559(subjIdx);
end



%% Compute residuals
additiveModelRes = additiveSensitivityRatioBestLambdaMax - realData(subjIdx);
IORes = IOSensitivityRatioBestLambdaMax - realData(subjIdx);

maxRes = max([abs(additiveModelRes); abs(IORes)]);

%% Compute RMSE
n = numel(realData(subjIdx));
AdditiveRMSE555 = sqrt((1/n).*sum((additiveSensitivityRatio555(subjIdx) - realData(subjIdx)).^2));
AdditiveRMSE559 = sqrt((1/n).*sum((additiveSensitivityRatio559(subjIdx) - realData(subjIdx)).^2));
AdditiveRMSE563 = sqrt((1/n).*sum((additiveSensitivityRatio563(subjIdx) - realData(subjIdx)).^2));

% compute additiveSensitivity ratio more finely

wavelengths = 550:0.5:570;

for i = 1:numel(wavelengths)
    [RedSensitivity, GreenSensitivity, ~] = simplePigmentModel(wavelengths(i));
    additiveSensitivityRatio = interp1(propLFine, RedSensitivity, propL)./interp1(propLFine, GreenSensitivity, propL);
    RMSE{subj}(i) = sqrt((1/n).*sum((additiveSensitivityRatio(subjIdx) - realData(subjIdx)).^2));

end



IORMSE555 = sqrt((1/n).*sum((IOSensitivityRatio555(subjIdx) - realData(subjIdx)).^2));
IORMSE559 = sqrt((1/n).*sum((IOSensitivityRatio559(subjIdx) - realData(subjIdx)).^2));
IORMSE563 = sqrt((1/n).*sum((IOSensitivityRatio563(subjIdx) - realData(subjIdx)).^2));





%% plot additive model vs real data

% 1:1 line
plot(ax1,[0 1.1*max(realData(:))],[0 1.1*max(realData(:))], 'LineWidth', markerLineWidth, 'LineStyle', '-', 'Color', defaultColor);
% plot data
plot(ax1, realData(subjIdx), additiveSensitivityRatioBestLambdaMax, 'LineStyle',...
    'none', 'Marker', AMmarker, 'MarkerSize', markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none', 'LineWidth', markerLineWidth);


%% plot additive model residuals

% % plot zero line
% plot(ax2, [0 1.05], [0 0], 'LineStyle', '-', 'LineWidth', markerLineWidth, 'Color', defaultColor)
% errorbar(ax2, propL(subjIdx), additiveModelRes, [], [], abs(minPropL(subjIdx)-propL(subjIdx)), abs(maxPropL(subjIdx)-propL(subjIdx)),...
%       'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color', markerEdgeColors(subj,:), 'capsize', 0);
% 
% plot(ax2, propL(subjIdx), additiveModelRes, 'LineStyle', 'none', 'Marker', AMmarker, ...
%     'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none', 'LineWidth', markerLineWidth);


%ylabel(additiveModelAx(subj), 'Residual');

%% plot IO vs real data
plot(ax1,realData(subjIdx),IOSensitivityRatioBestLambdaMax, 'LineStyle', 'none',...
    'Marker', IOmarker, 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none', 'LineWidth', markerLineWidth);

%% Plot IO residuals 
% 
% % plot zero line
% errorbar(ax2, propL(subjIdx), IORes, [], [], abs(minPropL(subjIdx)-propL(subjIdx)), abs(maxPropL(subjIdx)-propL(subjIdx)),...
%       'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color',  markerEdgeColors(subj,:), 'capsize', 0);
% plot(ax2,propL(subjIdx), IORes, 'LineStyle', 'none', 'Marker', ...
%     IOmarker, 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none', 'LineWidth', markerLineWidth);

%ylim(IOResAx(subj), [-1.1*maxRes 1.1*maxRes]);

%% plot io vs additive

%1:1
% plot(ax3,[0 1.1*max(realData(:))],[0 1.1*max(realData(:))], 'LineWidth', markerLineWidth, 'LineStyle', '-', 'Color', defaultColor);
% plot(ax3, additiveSensitivityRatioBestLambdaMax, IOSensitivityRatioBestLambdaMax,...
%     'LineStyle', 'none', 'Marker', defaultMarker, 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none', 'LineWidth', markerLineWidth);

%% RMSE for different lambda max

% additive model

RMSEaxLim = 1.1*max([AdditiveRMSE555; AdditiveRMSE559; AdditiveRMSE563; IORMSE555; IORMSE559; IORMSE563]);
plot(ax4, wavelengths, RMSE{subj}, 'LineStyle', '-', 'Color' ,markerEdgeColors(subj,:), 'LineWidth', markerLineWidth)
plot(ax4, [555.5 558.9 563.4], [AdditiveRMSE555 AdditiveRMSE559 AdditiveRMSE563], 'LineStyle', 'none',...
    'Color', markerEdgeColors(subj,:), 'LineWidth', markerLineWidth, 'Marker', AMmarker, 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none');


plot(ax4, [555.5 558.9 563.4], [IORMSE555 IORMSE559 IORMSE563], 'LineStyle', 'none',...
    'Color', markerEdgeColors(subj,:), 'LineWidth', markerLineWidth, 'Marker', IOmarker, 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColors(subj,:), 'MarkerFaceColor', 'none');
ylim(ax4, [0 RMSEaxLim]);


end