close all;
fontSize = 16;
markerLineWidth = 2;
axLineWidth = 2;
markerSize = 12;
markerEdgeColor = 'k';
markerFaceColor = [0.4 0.4 0.4; 0.9 0.9 0.9];
tickDir = 'in';

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
IO555 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO555.mat');
IO559 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO559.mat');
IO563 = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\IO563.mat');

IO_RedThresh555 = 10.^IO555.logThreshold680;
IO_GreenThresh555 = 10.^IO555.logThreshold543;
IO_RedThresh559 = 10.^IO559.logThreshold680;
IO_GreenThresh559 = 10.^IO559.logThreshold543;
IO_RedThresh563 = 10.^IO563.logThreshold680;
IO_GreenThresh563 = 10.^IO563.logThreshold543;

IOSensitivityRatio555 = IO_GreenThresh555./IO_RedThresh555;
IOSensitivityRatio559 = IO_GreenThresh559./IO_RedThresh559;
IOSensitivityRatio563 = IO_GreenThresh563./IO_RedThresh563;

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

RMSEfig = figure;
RMSEAxSolo = axes('Position', [14.75/20, 5.5/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth', axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'box', 'on'); hold(RMSEAxSolo, 'on');
grid(RMSEAxSolo, 'on');
ylim([0 0.015]);
xlabel(RMSEAxSolo, 'L-cone \lambda_{max}'); ylabel(RMSEAxSolo, 'RMSE');
set(RMSEfig, 'Units', 'inches')
set(RMSEfig, 'Position', [0.5   0.5    7    7])

set(RMSEAxSolo, 'Position', [0.175 0.175 0.8 0.8])

for subj = 1:numel(subjects)

% create figure for each subject
subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
subjFig(subj) = figure;

%% Set up axes

% Ideal observer vs real data axes
IOvsRealAx(subj) = axes('Position', [8/20, 12.25/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth',axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'yticklabels', '', 'box', 'on'); hold(IOvsRealAx(subj), 'on')
yticks(linspace(0, 0.03, 4));
xticks(linspace(0, 0.03, 4));
grid(IOvsRealAx(subj), 'on');
set(IOvsRealAx(subj), 'XTickLabelRotation', 0);
xlabel(IOvsRealAx(subj), 'S_{680}/S_{543} (obs.)');



% Ideal observer residual axes
IOResAx(subj) = axes('Position', [8/20, 5.5/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth',axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir,  'yticklabels', '','box', 'on'); hold(IOResAx(subj), 'on');
grid(IOResAx(subj), 'on');
xticks(0:(1/3):1);
xticklabels({'0', '0.33', '0.67', '1'});
set(IOResAx(subj),'XTickLabelRotation', 0);
ylim(IOResAx(subj), [-0.015 0.015]);
xlabel(IOResAx(subj), 'Proportion L'); 
xlim(IOResAx(subj), [0 1.05]); 

% Additve model vs real data axes
AdditiveVsRealAx(subj) =axes('Position', [3/20, 12.25/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth', axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'box', 'on'); hold(AdditiveVsRealAx(subj), 'on');
yticks(linspace(0, 0.03, 4));
xticks(linspace(0, 0.03, 4));
grid(AdditiveVsRealAx(subj), 'on');
set(AdditiveVsRealAx(subj), 'XTickLabelRotation', 0);
xlim(AdditiveVsRealAx(subj), [0 1.1*maxSensitivityRatio]); ylim(AdditiveVsRealAx(subj),[0 1.1*maxSensitivityRatio]);
xlim(IOvsRealAx(subj), [0 1.1*max(realData(:))]); ylim(IOvsRealAx(subj),[0 1.1*max(realData(:))]);
ylabel(AdditiveVsRealAx(subj), 'S_{680}/S_{543} (pred.)');
% Additve model residual axes
AdditiveResAx(subj) = axes('Position', [3/20, 5.5/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth', axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'box', 'on'); hold(AdditiveResAx(subj), 'on');
grid(AdditiveResAx(subj), 'on');
xticks(0:(1/3):1);
xticklabels({'0', '0.33', '0.67', '1'});
set(AdditiveResAx(subj),'XTickLabelRotation', 0);
xlim(AdditiveResAx(subj), [0 1.05]); 
%ylim(AdditiveResAx(subj), [-1.1*maxRes 1.1*maxRes]);
ylim(AdditiveResAx(subj), [-0.015 0.015]);
ylabel(AdditiveResAx(subj), 'Residual');

% Ideal observer vs additve model axes
IOvsAdditiveAx(subj) = axes('Position', [14.75/20, 12.25/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth', axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'box', 'on'); hold(IOvsAdditiveAx(subj), 'on');
grid(IOvsAdditiveAx(subj), 'on');
yticks(linspace(0, 0.03, 4));
xticks(linspace(0, 0.03, 4));
set(IOvsAdditiveAx(subj),'XTickLabelRotation', 0);
xlabel(IOvsAdditiveAx(subj), 'S_{680}/S_{543} (a.m.)');
ylabel(IOvsAdditiveAx(subj), 'S_{680}/S_{543} (i.o.)');
xlim(IOvsAdditiveAx(subj), [0 1.1*max(realData(:))]); ylim(IOvsAdditiveAx(subj),[0 1.1*max(realData(:))]);

% RMSE vs lambdaMax axes
RMSEAx(subj) = axes('Position', [14.75/20, 5.5/20, 4/20, 4/20], 'XColor', 'k', 'YColor', 'k', 'LineWidth', axLineWidth, 'FontSize', fontSize, 'TickDir', tickDir, 'box', 'on'); hold(RMSEAx(subj), 'on');
grid(RMSEAx(subj), 'on');
ylim([0 0.015]);
xlabel(RMSEAx(subj), 'L-cone \lambda_{max}'); ylabel(RMSEAx(subj), 'RMSE');
%% Get real data, IO, additive model predictions for subject

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
plot(AdditiveVsRealAx(subj),[0 1.1*max(realData(:))],[0 1.1*max(realData(:))], 'LineWidth', axLineWidth, 'LineStyle', '-', 'Color', 'k');
% plot data
plot(AdditiveVsRealAx(subj), realData(subjIdx), additiveSensitivityRatioBestLambdaMax, 'LineStyle', 'none', 'Marker', 'v', 'MarkerSize', markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:) , 'LineWidth', markerLineWidth);


%% plot additive model residuals

% plot zero line
plot(AdditiveResAx(subj), [0 1], [0 0], 'LineStyle', '-', 'LineWidth',axLineWidth, 'Color', 'k')
errorbar(AdditiveResAx(subj), propL(subjIdx), additiveModelRes, [], [], abs(minPropL(subjIdx)-propL(subjIdx)), abs(maxPropL(subjIdx)-propL(subjIdx)),...
      'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color', markerEdgeColor, 'capsize', 0);
  
plot(AdditiveResAx(subj), propL(subjIdx), additiveModelRes, 'LineStyle', 'none', 'Marker', 'v', ...
    'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:), 'LineWidth', markerLineWidth);


%ylabel(additiveModelAx(subj), 'Residual');

%% plot IO vs real data
% 1:1 line
plot(IOvsRealAx(subj),[0 1.1*max(realData(:))],[0 1.1*max(realData(:))], 'LineWidth', axLineWidth, 'LineStyle', '-', 'Color', 'k');
plot(IOvsRealAx(subj),realData(subjIdx),IOSensitivityRatioBestLambdaMax, 'LineStyle', 'none', 'Marker', '^', 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:), 'LineWidth', markerLineWidth);

%% Plot IO residuals 

% plot zero line
plot(IOResAx(subj), [0 1], [0 0], 'LineStyle', '-', 'LineWidth',axLineWidth, 'Color', 'k')
errorbar(IOResAx(subj), propL(subjIdx), IORes, [], [], abs(minPropL(subjIdx)-propL(subjIdx)), abs(maxPropL(subjIdx)-propL(subjIdx)),...
      'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color',  markerEdgeColor, 'capsize', 0);
plot(IOResAx(subj),propL(subjIdx), IORes, 'LineStyle', 'none', 'Marker', '^', 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:), 'LineWidth', markerLineWidth);

%ylim(IOResAx(subj), [-1.1*maxRes 1.1*maxRes]);

%% plot io vs additive

%1:1
plot(IOvsAdditiveAx(subj),[0 1.1*max(realData(:))],[0 1.1*max(realData(:))], 'LineWidth', axLineWidth, 'LineStyle', '-', 'Color', 'k');
plot(IOvsAdditiveAx(subj), additiveSensitivityRatioBestLambdaMax, IOSensitivityRatioBestLambdaMax,...
    'LineStyle', 'none', 'Marker', 'd', 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:), 'LineWidth', markerLineWidth);

%% RMSE for different lambda max

% additive model

RMSEaxLim = 1.1*max([AdditiveRMSE555; AdditiveRMSE559; AdditiveRMSE563; IORMSE555; IORMSE559; IORMSE563]);
plot(RMSEAx(subj), [555.5 558.9 563.4], [AdditiveRMSE555 AdditiveRMSE559 AdditiveRMSE563], 'LineStyle', 'none', 'Color', 'k', 'LineWidth', markerLineWidth, 'Marker', 'v', 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:));
plot(RMSEAx(subj), wavelengths, RMSE{subj}, 'LineStyle', '-', 'Color' ,'k', 'LineWidth', markerLineWidth)

plot(RMSEAx(subj), [555.5 558.9 563.4], [IORMSE555 IORMSE559 IORMSE563], 'LineStyle', 'none','Color', 'k', 'LineWidth', markerLineWidth, 'Marker', '^', 'MarkerSize',  markerSize, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:));
%ylim([0 RMSEaxLim]);
plot(RMSEAxSolo, [555.5 558.9 563.4], [AdditiveRMSE555 AdditiveRMSE559 AdditiveRMSE563], 'LineStyle', 'none', 'Color', 'k', 'LineWidth', markerLineWidth, 'Marker', 'v', 'MarkerSize', 18, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:));
plot(RMSEAxSolo, wavelengths, RMSE{subj}, 'LineStyle', '-', 'Color' ,'k', 'LineWidth', markerLineWidth)

plot(RMSEAxSolo, [555.5 558.9 563.4], [IORMSE555 IORMSE559 IORMSE563], 'LineStyle', 'none','Color', 'k', 'LineWidth', markerLineWidth, 'Marker', '^', 'MarkerSize',  18, 'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(subj,:));

set(subjFig(subj), 'Units', 'inches')
set(subjFig(subj), 'Position', [1 + subj   0    9    9])

end