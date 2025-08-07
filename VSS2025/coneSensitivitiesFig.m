% make cone fundamentals figure
axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
markerLineWidth = 1.5;
markerSize = 8;
tickLength = [0.025 0.025];

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors = FillInPhotoreceptors(photoreceptors);
wls = SToWls(photoreceptors.nomogram.S);

L = photoreceptors.quantalFundamentals(1,:);
M = photoreceptors.quantalFundamentals(2,:);
S = photoreceptors.quantalFundamentals(3,:);

figure; hold on

set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters', 'Color', 'w')
set(gca, 'Position', [0.2 0.2 0.8 0.8], 'FontSize', axFontSize, 'XColor', defaultColor, 'YColor', defaultColor, 'LineWidth', defaultLineWidth, 'TickDir', tickDir, 'box', 'off');
set(gca, 'TickLength', tickLength)

set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5   9   9])
set(gca, 'Position', [0.2 0.2 0.75 0.75])

plot(wls, log10(S), 'Color', [0 0 1], 'LineWidth',defaultLineWidth+1);
plot(wls, log10(M), 'Color', [0 1 0], 'LineWidth',defaultLineWidth+1);
plot(wls, log10(L), 'Color', [1 0 0], 'LineWidth',defaultLineWidth+1);

xlim([380 780]);
ylim([-6 0]);
xticks(300:100:800)
yticks(-6:2:0);

xlabel('Wavelength (nm)', 'FontSize', labelFontSize);
ylabel('Log sensitivity (normalized)', 'FontSize', labelFontSize);