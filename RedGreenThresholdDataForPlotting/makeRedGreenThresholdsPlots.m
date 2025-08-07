figure; hold on
gscatter(tbl.ProportionL, tbl.greenThresholdAUtimesGreenPower./tbl.redThresholdAUtimesRedPower, tbl.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Green threshold (au) * Green Power / Red threshold (au) * RedPower')
ylim([0 0.05])

figure; hold on
gscatter(tbl.ProportionL, tbl.incidentGreenQuanta./tbl.incidentRedQuanta, tbl.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Incident green quanta / incident red quanta')
ylim([0 0.05])

figure; hold on
gscatter(tbl.ProportionL, tbl.transmittedGreenQuanta./tbl.transmittedRedQuanta, tbl.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Transmitted green quanta / transmitted red quanta')
ylim([0 0.05])

