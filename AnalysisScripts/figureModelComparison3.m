close all;
clear;
axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
markers = {'^', 's', 'o'};
markerEdgeColors = [0 0 0.35; 1 0.4 0];
markerLineWidth = 1;
markerSize = 5;

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);
propL = ADT.L./(ADT.L + ADT.M);
minPropL = ADT.L./(ADT.L + ADT.M + ADT.X);
maxPropL = (ADT.L + ADT.X)./(ADT.L + ADT.M + ADT.X);

convertThresholdsToQuanta;

figure;
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5    13.5   13.5], 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters');


ax1 = axes('Position', [0.125 0.4, 0.35 0.35], 'LineWidth', defaultLineWidth, 'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
hold(ax1, 'on');
xlabel(ax1, 'Obs. sensitivity ratio', 'FontSize', labelFontSize);
ylabel(ax1, 'Pred. sensitivity ratio', 'FontSize', labelFontSize);
xlim(ax1, [0 0.035]); ylim(ax1, [0 0.035]);

realDataMax = max(retIrradGreenQuantaPerSecM2./retIrradRedQuantaPerSecM2);
plot(ax1,[0 1.1*realDataMax],[0 1.1*realDataMax], 'LineWidth',...
    markerLineWidth, 'LineStyle', '-', 'Color', defaultColor);

ax2 = axes('Position', [0.625, 0.4, 0.35 0.35], 'LineWidth', defaultLineWidth,  'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', axFontSize, 'TickDir',tickDir);
hold(ax2, 'on')
xlabel(ax2, '\lambda_{max} (nm)', 'FontSize', labelFontSize);
ylabel(ax2, 'RMSE', 'FontSize', labelFontSize);
xlim(ax2, [550 570]);

coneSpectralPeaks = [555.5 530.3 420.7; 558.9 530.3 420.7; 563.4 530.3 420.7];

% for each subject
% for each cone spectral peak

% get model data at tested propL and at wavelengths of interest

for subj = 1:numel(subjects)
    subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
    params = getSubjectParams(subjects{subj});
    subjPropL = propL(subjIdx);
    realData{subj} = retIrradGreenQuantaPerSecM2(subjIdx)./retIrradRedQuantaPerSecM2(subjIdx);
    n = numel(realData{subj});
    for p = 1:size(coneSpectralPeaks,1)
        params.nomogram.lambdaMax = transpose(coneSpectralPeaks(p,:));
        [RedSensitivity, GreenSensitivity] = simplePigmentModel(params, subjPropL);
        additiveSensitivityRatio{subj,p} = RedSensitivity./GreenSensitivity;
        additiveRMSE{subj,p} = sqrt((1/n).*sum((additiveSensitivityRatio{subj,p} - realData{subj}).^2));

        if subj == 1 && p == 3 || subj == 2 && p == 2
            plot(ax1, realData{subj}, additiveSensitivityRatio{subj,p},...
                'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'o', 'MarkerSize', markerSize, 'MarkerEdgeColor', ...
                markerEdgeColors(subj,:), 'MarkerFaceColor', markerEdgeColors(subj,:));
        else
          %  plot(ax1, realData{subj}, additiveSensitivityRatio{subj,p}, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', markerEdgeColors(subj,:))
        end
    end
end

wavelengths = 550:0.5:570;
RMSEmax = 0;
for subj = 1:numel(subjects)
    subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
    params = getSubjectParams(subjects{subj});
    subjPropL = propL(subjIdx);
    realData{subj} = retIrradGreenQuantaPerSecM2(subjIdx)./retIrradRedQuantaPerSecM2(subjIdx);
    n = numel(realData{subj});
    for wl = 1:numel(wavelengths)
        params.nomogram.lambdaMax = [wavelengths(wl) 530.3 420.7]';
        [RedSensitivity, GreenSensitivity] = simplePigmentModel(params, subjPropL);
        additiveSensitivityRatioFine = RedSensitivity./GreenSensitivity;
        additiveRMSEFine{subj,wl} = sqrt((1/n).*sum((additiveSensitivityRatioFine - realData{subj}).^2));

        RMSEmax = max([RMSEmax additiveRMSEFine{subj,wl}]);
    end
    plot(ax2, wavelengths, [additiveRMSEFine{subj,1:end}], 'LineWidth', markerLineWidth, 'Color', markerEdgeColors(subj,:))
    plot(coneSpectralPeaks(:,1), [additiveRMSE{subj, 1:end}], 'LineStyle', ...
        'none', 'LineWidth', markerLineWidth, 'Marker', 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', markerEdgeColors(subj,:), 'MarkerEdgeColor', markerEdgeColors(subj,:))
end
ylim([0 RMSEmax]);
ax2.YAxis.Exponent = -2;