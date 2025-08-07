% for axes
fontSize = 12;
lineWidth = 2;
markerSize = 3;
defaultColor = [0 0 0];
tickDir = 'out';
tickLength = [0.025 0.025];

% load cone apertures
CA_10001R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\OCTMaster\coneApertures.mat');
CA_20217R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\20217R\OCTMaster\coneApertures.mat');
CA_10001R = CA_10001R./max(CA_10001R(:));
CA_20217R = CA_20217R./max(CA_20217R(:));

OCT_10001R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\OCTMaster\Raw_image_fundus_scaled.tif');
OCT_20217R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\20217R\OCTMaster\Raw_image_fundus_scaled.tif');
coneData10001R = load('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\OCTMaster\coneData.mat', 'coneData');
coneData20217R = load('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\20217R\OCTMaster\coneData.mat', 'coneData');

coneData10001R = coneData10001R.coneData;
coneData20217R = coneData20217R.coneData;

likelyMissingCones10001R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\OCTMaster\likelyMissingCones.mat');
likelyMissingCones20217R = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\20217R\OCTMaster\likelyMissingCones.mat');

% Get targeted locations for each subject

[xm_10001R, ym_10001R] = findTargetedLocationsInMaster('10001R');
[xm_20217R, ym_20217R] = findTargetedLocationsInMaster('20217R');

% find horizontal and vertical extents of stimulated regions in each
% subject

xRange10001R = max(xm_10001R) - min(xm_10001R);
yRange10001R = max(ym_10001R) - min(ym_10001R);
xRange20217R = max(xm_20217R) - min(xm_20217R);
yRange20217R = max(ym_20217R) - min(ym_20217R);

xRange = max([xRange10001R xRange20217R]);
yRange = max([yRange10001R yRange20217R]);

rectHeight = max([xRange yRange]) + 42;
rectWidth = max([xRange yRange]) +42;


rect10001R = [220 239 rectWidth+1 rectHeight+1];
rect20217R = [min(xm_20217R)-21, min(ym_20217R)-21, rectWidth, rectHeight+1];

% CA_10001R_cropped = imcrop(CA_10001R, rect10001R);
% CA_20217R_cropped = imcrop(CA_20217R, rect20217R);
OCT_10001R_cropped = imcrop(OCT_10001R, rect10001R);
OCT_20217R_cropped = imcrop(OCT_20217R, rect20217R);

cone_locs_10001R_cropped = str2double(coneData10001R(:, 1:2)) - [220-1 239-1];
cone_locs_20217R_cropped = str2double(coneData20217R(:, 1:2)) - [(min(xm_20217R)- 21-1) (min(ym_20217R)-21)-1];

badIdx10001R = any(cone_locs_10001R_cropped <= 0, 2) | cone_locs_10001R_cropped(:,1) > size(OCT_10001R_cropped,2) | cone_locs_10001R_cropped(:,2) > size(OCT_10001R_cropped,1);
badIdx20217R = any(cone_locs_20217R_cropped <= 0, 2) | cone_locs_20217R_cropped(:,1) > size(OCT_20217R_cropped,2) | cone_locs_20217R_cropped(:,2) > size(OCT_20217R_cropped,1);

xm_10001R_cropped = xm_10001R - 220 - 1;
ym_10001R_cropped = ym_10001R - 239 - 1;

xm_20217R_cropped = xm_20217R - min(xm_20217R)+21 - 1;
ym_20217R_cropped = ym_20217R- min(ym_20217R)+21 - 1;

likelyMissingCones10001R_cropped = likelyMissingCones10001R - [220-1 239-1];
likelyMissingCones20217R_cropped = likelyMissingCones20217R - [(min(xm_20217R)- 21-1) (min(ym_20217R)-21)-1];


% Now make figure

h = figure;
set(h, 'Units', 'centimeters');
set(h, 'Position', [1 1 13.5 13.5]);
set(h, 'Renderer', 'painters', 'PaperPositionMode', 'auto', 'Color', 'w');


axMosaic10001R = axes('Position', [0.1 0.5 0.425 0.425], 'Ydir', 'reverse');xlim([1 size(OCT_10001R_cropped,2)]), ylim([1 size(OCT_10001R_cropped,1)]); hold on; axis off;
axMosaic20217R = axes('Position', [0.55 0.5 0.425 0.425], 'ydir', 'reverse');xlim([1 size(OCT_20217R_cropped,2)]), ylim([1, size(OCT_10001R_cropped,1)]); hold on; axis off

imagesc(axMosaic10001R, OCT_10001R_cropped); colormap gray, 
imagesc(axMosaic20217R, OCT_20217R_cropped); colormap gray, 

% plot cone locs
LColor = [1 0 0];
MColor = [0 0.8 0];
SColor = [0 0 1];

plot(axMosaic10001R, cone_locs_10001R_cropped(coneData10001R(:,3) == "L" & ~badIdx10001R,1), cone_locs_10001R_cropped(coneData10001R(:,3) == "L" & ~badIdx10001R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', LColor, 'MarkerSize', markerSize)
plot(axMosaic10001R,cone_locs_10001R_cropped(coneData10001R(:,3) == "M" & ~badIdx10001R,1), cone_locs_10001R_cropped(coneData10001R(:,3) == "M" & ~badIdx10001R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', MColor, 'MarkerSize', markerSize)
plot(axMosaic10001R, cone_locs_10001R_cropped(coneData10001R(:,3) == "S" & ~badIdx10001R,1), cone_locs_10001R_cropped(coneData10001R(:,3) == "S" & ~badIdx10001R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', SColor, 'MarkerSize', markerSize)
% not L, M, or S index (i.e. missing)
labeledMissing10001R = coneData10001R(:,3) ~= "L" & coneData10001R(:,3) ~= "M" & coneData10001R(:,3) ~= "S";
%plot(axMosaic10001R, cone_locs_10001R_cropped(labeledMissing10001R & ~badIdx10001R,1), cone_locs_10001R_cropped(labeledMissing10001R & ~ badIdx10001R,2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'y', 'MarkerSize', 3);



plot(axMosaic20217R, cone_locs_20217R_cropped(coneData20217R(:,3) == "L" & ~badIdx20217R,1), cone_locs_20217R_cropped(coneData20217R(:,3) == "L" & ~badIdx20217R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', LColor, 'MarkerSize', markerSize)
plot(axMosaic20217R,cone_locs_20217R_cropped(coneData20217R(:,3) == "M" & ~badIdx20217R,1), cone_locs_20217R_cropped(coneData20217R(:,3) == "M" & ~badIdx20217R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', MColor, 'MarkerSize', markerSize)
plot(axMosaic20217R, cone_locs_20217R_cropped(coneData20217R(:,3) == "S" & ~badIdx20217R,1), cone_locs_20217R_cropped(coneData20217R(:,3) == "S" & ~badIdx20217R,2),...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', SColor, 'MarkerSize', markerSize)
labeledMissing20217R = coneData20217R(:,3) ~= "L" & coneData20217R(:,3) ~= "M" & coneData20217R(:,3) ~= "S";
%plot(axMosaic20217R, cone_locs_20217R_cropped(labeledMissing20217R & ~badIdx20217R,1), cone_locs_20217R_cropped(labeledMissing20217R & ~ badIdx20217R,2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'y', 'MarkerSize', 3);

missing10001R = [cone_locs_10001R_cropped(labeledMissing10001R,:) ;likelyMissingCones10001R_cropped];
missing20217R = [cone_locs_20217R_cropped(labeledMissing20217R,:) ;likelyMissingCones20217R_cropped];

% Plot likely missing cones if they are within stimulated region

 for i = 1:size(missing10001R,1)
    %x = likelyMissingCones10001R_cropped(i,1); y = likelyMissingCones10001R_cropped(i,2);
    if any(all(abs(missing10001R(i,:) - [xm_10001R_cropped' ym_10001R_cropped']) <= 11, 2)) &  all(any(abs(missing10001R(i,:) - cone_locs_10001R_cropped) > 2, 2))
        plot(axMosaic10001R, missing10001R(i,1), missing10001R(i,2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'y', 'MarkerSize', markerSize)
    else
    end


end

for i = 1:size(missing20217R,1)
    %x = likelyMissingCones10001R_cropped(i,1); y = likelyMissingCones10001R_cropped(i,2);
    if any(all(abs(missing20217R(i,:) - [xm_20217R_cropped' ym_20217R_cropped']) <= 11, 2)) &  all(any(abs(missing20217R(i,:) - cone_locs_20217R_cropped) > 2, 2))
        plot(axMosaic20217R, missing20217R(i,1), missing20217R(i,2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'y', 'MarkerSize', markerSize)
    else
    end


end



for i = 1:numel(xm_10001R)
    drawrectangle(axMosaic10001R, 'Position', [xm_10001R_cropped(i)-11 ym_10001R_cropped(i)-11 21 21], 'LineWidth', 1,  'Color', 'w',  'FaceAlpha', 0.2, 'InteractionsAllowed', 'none');
end

for i = 1:numel(xm_20217R)
    drawrectangle(axMosaic20217R, 'Position', [xm_20217R_cropped(i)-11 ym_20217R_cropped(i)-11 21 21], 'LineWidth', 1, 'Color', 'w' , 'FaceAlpha', 0.2, 'InteractionsAllowed', 'none');
end


AssessRandomnessOfObservedPropL;

axHist10001R = copyobj(histAxes(1), h);
axHist10001R.Position = [0.1 0.325 0.425 0.15];
axHist20217R = copyobj(histAxes(2), h);
axHist20217R.Position = [0.55 0.325 0.425 0.15];
title(axHist10001R, ''); title(axHist20217R, '');
ylabel(axHist20217R, '');
axHist20217R.YTickLabel = '';


axMosaic10001R.Position = [0.1 0.5 0.425 0.425];
axMosaic20217R.Position = [0.55 0.5 0.425 0.425];

% make axes look nice
set(axHist10001R, 'XColor', defaultColor, 'YColor', defaultColor, 'FontSize', fontSize, 'LineWidth', lineWidth, 'TickDir', tickDir, 'TickLength', tickLength)
set(axHist20217R, 'XColor', defaultColor, 'YColor', defaultColor,'FontSize', fontSize, 'LineWidth', lineWidth, 'TickDir', tickDir, 'TickLength', tickLength)
xlim(axHist10001R, [0 1]);
xlim(axHist20217R, [0 1]);

% 
% [RedSensitivity555, GreenSensitivity555, propLFine] = simplePigmentModel(555.5);
% [RedSensitivity559, GreenSensitivity559, ~] = simplePigmentModel;
% [RedSenstivity563, GreenSensitivity563, ~] = simplePigmentModel(563.4);
% 
% axSimpleModel = axes(h,'Position', [13.75/20 7.5/20 6/20 7/20], 'YColor', 'k', 'XColor', 'k', 'FontSize', fontSize, 'LineWidth', lineWidth); hold(axSimpleModel, 'on');
% ylabel(axSimpleModel, 'S_{680}/S_{543}');
%ylim(axSimpleModel, [0 0.035]);

% 
% plot(axSimpleModel, propLFine, RedSensitivity555./GreenSensitivity555, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':'); 
% plot(axSimpleModel, propLFine, RedSensitivity559./GreenSensitivity559, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-');
% plot(axSimpleModel, propLFine, RedSenstivity563./GreenSensitivity563, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');

