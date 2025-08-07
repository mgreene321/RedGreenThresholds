% make demo cone mosaic image

fontSize = 12;
lineWidth = 2;
markerSize = 4.5;
defaultColor = [0 0 0];

cone_data = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\5_21_2024\TransferFromOCT\dailyData.mat');
cone_locs = str2double(dailyData(:, [1 2])); % N x 2 matrix of cone coordinates
cone_labels = dailyData(:,3);


filePath = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\ConeTransfersTo5_21_2024';
retina_map = importdata(fullfile(filePath, 'meanMeanFrame.mat'));

figure, imshow(retina_map, [min(retina_map(:)) max(retina_map(:))]); hold on
set(gcf, 'Units', 'centimeters');
%set(gcf, 'Position', [1 1 18 18]);
set(gcf, 'Renderer', 'painters', 'PaperPositionMode', 'auto', 'Color', 'w');
LColor = [1 0 0];
MColor = [0 0.8 0];
SColor = [0 0 1];

plot(cone_locs(cone_labels=='L',1), cone_locs(cone_labels=='L',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', LColor, 'MarkerSize', markerSize); 
plot(cone_locs(cone_labels=='M',1), cone_locs(cone_labels=='M',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', MColor, 'MarkerSize', markerSize); 
plot(cone_locs(cone_labels=='S',1), cone_locs(cone_labels=='S',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', SColor, 'MarkerSize', markerSize); 