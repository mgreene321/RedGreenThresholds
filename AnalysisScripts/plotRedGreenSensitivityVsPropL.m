convertThresholdsToQuanta

greenFWHM = 26.2871; % nm
redFWHM = 30.6503;
greenSigma = greenFWHM/(2*sqrt(2*log(2)));
redSigma = redFWHM/(2*sqrt(2*log(2)));
Red = normpdf(wvl, 680, redSigma);
Green = normpdf(wvl,543, greenSigma);
Red = Red./max(Red(:));
Green = Green./max(Green(:));
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
propL = ADT.L./(ADT.L + ADT.M); 
fontSize = 14;
%% Plot red:green sensitivity vs propL
figure; hold on; set(gca, 'FontSize', fontSize);


xlabel('Proportion L')
ylabel('Red : Green Sensitivity');
markerFaceColors = {'k', 'w'};

[RedSensitivity555, GreenSensitivity555, propLFine] = simplePigmentModel(555.5);
[RedSensitivity559, GreenSensitivity559, ~] = simplePigmentModel;
[RedSenstivity563, GreenSensitivity563, ~] = simplePigmentModel(563.4);

%axSimpleModel = axes(h,'Position', [12.5/20 7.5/20 7/20 7/20]); hold(axSimpleModel, 'on');
plot(propLFine, RedSensitivity555./GreenSensitivity555, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':'); 
plot(propLFine, RedSensitivity559./GreenSensitivity559, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
plot(propLFine, RedSenstivity563./GreenSensitivity563, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');


for s = 1:length(subjects)
subjectId = subjects{s};
subjectIdx = strcmpi(ADT.SubjectID, subjectId); 

plot(propL(subjectIdx), incidentGreenQuanta(subjectIdx)./incidentRedQuanta(subjectIdx), 'LineStyle', 'none', 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', markerFaceColors{s}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10)

end
% 
% RedSensitivity = propLFine .* (dot(Red', L) - dot(Red', M)) + dot(Red', M);
% GreenSensitivity = propLFine .* (dot(Green', L) - dot(Green', M)) + dot(Green', M);
% plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% legend show
