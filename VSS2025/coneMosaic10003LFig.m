% make demo cone mosaic image

fontSize = 12;
lineWidth = 2;
markerSize = 2.5;
defaultColor = [0 0 0];


fileName ='output_04;19;23-19;31_optimal_2k2k_region.json'; %'output_01;12;23-15;07_optimal_2k2k_region.json';
%'
filePath = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\Unstabilized1deg\10003L\size30dur1_10_27_2022_11_4_51\ConeTransfer';
%'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\Unstabilized1deg\20236R\size100dur1_9_20_2022_15_10_29\ConeTransfer';

str = fileread(fullfile(filePath, fileName));
data = jsondecode(str);
cone_data = transpose([data.cone_data.l]);
cone_locs = cone_data(:, [1 2]); %N x 2 matrix, rows are cones, columns x and y coords
cone_labels = [data.cone_data.c]'; %N x 1 matrix with cone types
%retina_map = imread('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\Unstabilized1deg\20236R\size100dur1_9_20_2022_15_10_29\retina_map.png');
 
%retina_map =  imread('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\Unstabilized1deg\10003L\size30dur1_10_27_2022_11_4_51\retina_map.png');
load('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\Unstabilized1deg\10003L\size30dur1_10_27_2022_11_4_51\rslamv2_outputs.mat');
figure, imshow(retina_map, [min(retina_map(:)) max(retina_map(:))]); hold on
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 1 18 18]);
set(gcf, 'Renderer', 'painters', 'PaperPositionMode', 'auto', 'Color', 'w');
LColor = [1 0 0];
MColor = [0 0.8 0];
SColor = [0 0 1];

plot(cone_locs(cone_labels=='l',1), cone_locs(cone_labels=='l',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', LColor, 'MarkerSize', markerSize); 
plot(cone_locs(cone_labels=='m',1), cone_locs(cone_labels=='m',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', MColor, 'MarkerSize', markerSize); 
plot(cone_locs(cone_labels=='s',1), cone_locs(cone_labels=='s',2), 'LineStyle', ...
    'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', SColor, 'MarkerSize', markerSize); 