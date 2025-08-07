close all; clear;
axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
markerLineWidth = 1.5;
markerSize = 8;
tickLength = [0.025 0.025];
lineStyles = {':', '-', '--'};

convertThresholdsToQuanta

ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
propL = ADT.L./(ADT.L + ADT.M); 
minPropL = ADT.L./(ADT.L + ADT.M + ADT.X);
maxPropL = (ADT.L + ADT.X)./(ADT.L + ADT.M + ADT.X);
subjects = unique(ADT.SubjectID);

%% Plot red:green sensitivity vs propL
figure; hold on;
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters', 'Color', 'w')
set(gca, 'Position', [0.225 0.2 0.8 0.8], 'FontSize', axFontSize, 'XColor', defaultColor, 'YColor', defaultColor, 'LineWidth', defaultLineWidth, 'TickDir', tickDir, 'box', 'off');
set(gca, 'TickLength', tickLength)

xticks(0:(1/3):1), xticklabels({'0', '0.33', '0.67', '1'});
yticks(linspace(0, 0.03, 4));

xlabel('Proportion L', 'FontSize', labelFontSize)
ylabel('Sensitivity ratio (680/543)', 'FontSize', labelFontSize);
markerFaceColors = [0 0 0.35; 1 0.4 0];

coneSpectralPeaks = [555.5 530.3 420.7; 558.9 530.3 420.7; 563.4 530.3 420.7];
 propLFine = linspace(0,1,1e3);
  % optimal L cone weight
        optimal_wL = findOptimalConeWeights;
for subj = 1:numel(subjects)
    params = getSubjectParams(subjects{subj});
    for p = 1:size(coneSpectralPeaks,1)
        params.nomogram.lambdaMax = transpose(coneSpectralPeaks(p,:));
       
       

        [RedSensitivity{subj,p}, GreenSensitivity{subj,p}] = simplePigmentModel(params, propLFine, optimal_wL{subj,p});
        plot(propLFine, RedSensitivity{subj,p}./GreenSensitivity{subj,p}, 'LineWidth', 1.5, 'Color', markerFaceColors(subj,:), 'LineStyle', lineStyles{p}, 'HandleVisibility', 'off');

    end
end


%axSimpleModel = axes(h,'Position', [12.5/20 7.5/20 7/20 7/20]); hold(axSimpleModel, 'on');



% 

for i = 1:size(ADT,1)
     subjectId = ADT.SubjectID{i};
     subjectIdx = strcmpi(subjects, subjectId); 
    errorbar(propL(i), retIrradGreenQuantaPerSecM2(i)./retIrradRedQuantaPerSecM2(i), [], [], abs(minPropL(i)-propL(i)), abs(maxPropL(i)-propL(i)), ...
        'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color',defaultColor, 'capsize', 0);
    
end
for s = 1:length(subjects)
subjectId = subjects{s};
subjectIdx = strcmpi(ADT.SubjectID, subjectId); 

errorbar(propL(subjectIdx), 10.^log10(retIrradGreenQuantaPerSecM2(subjectIdx)./retIrradRedQuantaPerSecM2(subjectIdx)), [], [], abs(minPropL(subjectIdx)-propL(subjectIdx)), abs(maxPropL(subjectIdx)-propL(subjectIdx)), ...
        'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', 'none', 'Color', defaultColor, 'capsize', 0);


plot(propL(subjectIdx), 10.^log10(retIrradGreenQuantaPerSecM2(subjectIdx)./retIrradRedQuantaPerSecM2(subjectIdx)),...
    'LineStyle', 'none', 'LineWidth', markerLineWidth, 'Marker', defaultMarker, 'MarkerFaceColor', markerFaceColors(s,:), 'MarkerEdgeColor', defaultColor, 'MarkerSize', markerSize, 'DisplayName', subjects{s})

end


xlim([0 1.05]);
ylim([0 1.1*max(retIrradGreenQuantaPerSecM2./retIrradRedQuantaPerSecM2)]);
%xticks([0 1/3 2/3 1])
%xticklabels({'0', '0.33', '0.67', '1'})

set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5    9    9])

set(gca, 'Position', [0.2 0.2 0.75 0.75])
% 
% RedSensitivity = propLFine .* (dot(Red', L) - dot(Red', M)) + dot(Red', M);
% GreenSensitivity = propLFine .* (dot(Green', L) - dot(Green', M)) + dot(Green', M);
% plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% legend show
