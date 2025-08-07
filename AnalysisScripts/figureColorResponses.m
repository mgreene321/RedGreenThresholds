axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
markerLineWidth = 1.5;
tickLength =[0.025 0.025];


%% set up axes
figure; set(gcf, 'renderer', 'painters', 'paperpositionmode', 'auto', 'color', 'w');
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [1 1 18 18]);

axTopRight = axes('Position', [0.7125 0.6 0.25 0.25]);
axTopMid = axes('Position', [0.4375 0.6 0.25 0.25]);
axBottomRight = axes('Position', [0.7125 0.325 0.25 0.25]);
axBottomMid = axes('Position', [0.4375 0.325 0.25 0.25]);
axTopLeft = axes('Position', [0.125 0.6 0.1875 0.25]);
axBottomLeft = axes('Position', [0.125 0.325 0.1875 0.25]);

colorMatAxes = {axTopRight, axTopMid, axBottomRight, axBottomMid};
barAxes = {axTopLeft, axBottomLeft};
for i = 1:numel(colorMatAxes)
    hold(colorMatAxes{i}, 'on')
    colorMatAxes{i}.FontSize = axFontSize;
    colorMatAxes{i}.XColor = defaultColor; colorMatAxes{i}.YColor = defaultColor;
    colorMatAxes{i}.TickDir = tickDir;
    colorMatAxes{i}.Box = 'on';
    colorMatAxes{i}.LineWidth = defaultLineWidth;
    colorMatAxes{i}.GridLineWidth = 1.25;
end



% Load in data table
%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
% Eye constants

pupilDiamMm = 6.5;
eyeLengthMm = 16.7;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);

stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
stimSideLengthPix = 21;
pixPerDeg = 560;
stimSideLengthDeg = stimSideLengthPix/pixPerDeg;
stimSideLengthUm = DegreesToRetinalMM(eyeLengthMm, stimSideLengthDeg) * 1e3;
stimAreaUm2 = stimSideLengthUm^2;
%ISareaUm2 = pi*(ISdiameter/2)^2;

coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin


root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);
%RDT = RDT(strcmpi(RDT.SubjectID, subjectId),:);
% Load in aanalyzed table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(RDT.SubjectID);

% convert from intensity in au to "transmitted quanta" (quantal flux at
% photoreceptors)


% group table by folder
photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors = FillInPhotoreceptors(photoreceptors);
scanAngleDeg = 0.9;

%% Account for bleaching by projector background


indices = getLocationIndicesInRawDataTable;
propThresh = nan(size(RDT,1),1);
for i = 1:size(RDT,1)
    temp = cellfun(@(x) find(x), indices, 'UniformOutput', false);
    ADT_idx = cellfun(@(x) ismember(i, x), temp);
    if RDT.Channel(i) == 1

        propThresh(i,:) = RDT.IntensityAU(i)./ADT.RedThresholdPowerAU(ADT_idx);

    elseif RDT.Channel(i) == 2

        propThresh(i,:) = RDT.IntensityAU(i)./ADT.GreenThresholdPowerAU(ADT_idx);
    end

end
intensity = propThresh;
propLEdges = linspace(0,1, 9);


minIntensity= min(intensity(RDT.YesNo==1));
maxIntensity= max(intensity(RDT.YesNo==1));

%intEdges = 10.^(log10(minIntensity):log10(2):log10(maxIntensity));
intEdges = 2.^(-2:5);

%% colormaps
numColorMapRows = 256;

rgbMatrix_Red = [linspace(0,1, 256)', linspace(0, 0.25, 256)', linspace(0,0.25, 256)'];
rgbMatrix_Green = [linspace(0, 0.33, 256)', linspace(0,0.9, 256)', linspace(0, 0.33, 256)'];

colormaps = cat(3, rgbMatrix_Red, rgbMatrix_Green);

propL = RDT.L./(RDT.L + RDT.M);
minTrialsPerBin =10;
colorMatCounter = 1;


for subj = 1:numel(subjects)
    subjectId = subjects{subj};
    for ch = 1:2
        for col = 1:4

            idx = strcmpi(RDT.SubjectID, subjectId) & RDT.Channel == ch & RDT.Color == col;
            N{ch, col, subj} = histcounts2(propL(idx), intensity(idx),propLEdges, intEdges);

        end
    end
end


for subj = 1:numel(subjects)
    for ch = 1:2
        Nseen{ch, subj} = sum(cat(3, N{ch, 1:3, subj}),3);
        %Nseen{ch, subj} = Nseen{ch, subj} .* (Nseen{ch,subj} >= minTrialsPerBin);
        NexpectedColor{ch, subj} = N{ch,ch, subj};
        NunexpectedColor{ch, subj} = sum(cat(3, N{ch, 1:2, subj}),3) - NexpectedColor{ch,subj};
        temp = NexpectedColor{ch, subj}./Nseen{ch, subj};
        % handle nan and inf

        temp = log2(temp);
        temp(temp == -Inf) = -10;
        temp(~isfinite(temp)) = nan;
        temp = rot90(temp);
        temp = flipud(temp);
        colorMat{ch, subj} = temp;

        %unexpected color
        temp2 = NunexpectedColor{ch,subj}./sum(cat(3, N{ch, 1:3, subj}),3);
        %temp2 = log2(temp2);
        temp2(temp2 == -Inf) = -10;
        temp2(~isfinite(temp2)) = nan;
        temp2 = rot90(temp2);
        temp2 = flipud(temp2);
        colorMat2{ch,subj} = temp2;



    end
end


for subj = 1:numel(subjects)

    for ch = 1:2
        currentAxes = colorMatAxes{colorMatCounter};

        numRows = size(colorMat{ch, subj},1);
        numCols = size(colorMat{ch, subj},2);
        xVals = 0:1:numCols;
        yVals = 0:1:numRows;
        [xGrid, yGrid] = meshgrid(xVals, yVals);
        linesPerSquareEdge = 2;
        numBoxes = max([numCols numRows]); % Make it square for the largest dimension of the data matrix

       x1 = [zeros(1,linesPerSquareEdge*(numBoxes+1)); linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))];
       y1 =  [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))];
        x2 = [zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes; linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))];
        y2 =  [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes];

        x1(:,end) = [];
        y1(:,end) = [];
        x2(:,end) = [];
        y2(:,end) = [];

        xx = [x1 x2];
        yy = [y1 y2];
        plot(currentAxes, xx, yy, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 1);

      

        % plot(currentAxes, [zeros(1,linesPerSquareEdge*(numBoxes+1)); linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))],...
        %     [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))], 'k-'); hold on;
        % plot(currentAxes, [zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes; linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))],...
        %     [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes], 'k-'); hold on;

        for r = 1:numRows
            for c = 1:numCols
                if isfinite(colorMat{ch,subj}(r,c))
                    patch(currentAxes, [xGrid(r,c:c+1) fliplr(xGrid(r+1,c:c+1))], [yGrid(r,c:c+1) fliplr(yGrid(r+1,c:c+1))], colorMat{ch, subj}(r,c), 'EdgeColor', 'w', 'linewidth', 1.25);
                 
                end
            end
        end

        xlim(currentAxes, [0 numel(propLEdges)-1]);
        ylim(currentAxes, [0 numel(intEdges)-1]);

        propLTicks = 0+mod(numel(propLEdges)-1, 2):2: numel(propLEdges) -mod(numel(propLEdges)-1, 2);
        intTicks = 0+mod(numel(intEdges)-1, 2):2:numel(intEdges) -mod(numel(intEdges)-1, 2);

        % xticks(currentAxes, propLTicks);
        % yticks(currentAxes, intTicks);
        xticks(currentAxes, 0:1:numel(propLEdges));
        yticks(currentAxes, 0:1:numel(intEdges));

         grid(currentAxes, 'on');
         set(currentAxes, 'GridColor', defaultColor, 'GridAlpha', 1, 'GridLineStyle', ':')

        xTickLabs =  cellfun(@(X) num2str(X, '%0.3g'), num2cell(propLEdges), 'UniformOutput', false);
        tempXTickLabs = cell(size(xTickLabs));
        tempXTickLabs(:) = {''};
        tempXTickLabs(1:4:end) = xTickLabs(1:4:end);
        xTickLabs = tempXTickLabs;


        yTickLabs =  cellfun(@(X) num2str(X, '%0.3g'), num2cell(round(intEdges,2)), 'UniformOutput', false);

        xTickLabs(2:2:end) = {''};
        yTickLabs(1:2:end) = {''};

        xticklabels(currentAxes, xTickLabs);
        yticklabels(currentAxes, yTickLabs);
        currentAxes.XTickLabelRotation = 0;
        
       


        currentAxes.Colormap = colormaps(:,:,ch);
clim(currentAxes, [log2(0.0625/1.25) log2(1)]);
        if colorMatCounter <=2

            currentAxes.XTickLabel  = '';
        cb{colorMatCounter} = colorbar(currentAxes, 'LineWidth', defaultLineWidth, 'FontSize', axFontSize, 'TickDirection', tickDir, 'Color', defaultColor, 'Location', 'northoutside');
        %caxis([0 1.01]);

        PexpectedTicks= 2.^(log2(0.0625):log2(2):log2(1));
        %   clim(cb{colorMatCounter},[log2(0.0625/1.25) log2(1.25)]);
        cb{colorMatCounter}.Limits = [log2(0.0625/1.25) log2(1)];
        cb{colorMatCounter}.Ticks = [-4:1:0];
        cb{colorMatCounter}.TickLabels = cellfun(@(X) [num2str(X, '%0.3g') '%'], num2cell(100.*PexpectedTicks), 'UniformOutput', false);
        else
            xlabel(currentAxes, 'Proportion L', 'FontSize', labelFontSize);
        end

        if mod(colorMatCounter,2) ~= 0
            currentAxes.YTickLabel = '';
        else
             ylabel(currentAxes, 'Intensity (x threshold)', 'FontSize', labelFontSize);
        end


        
        colorMatCounter = colorMatCounter+1;

    end
end

for i = 1:numel(colorMatAxes)
    colorMatAxes{i}.Position(4) = 0.25;
end

cb{2}.Position(2) = cb{1}.Position(2);

%% Histograms

xpositions = [1/3 2/3];
for s = 1:numel(subjects)
    for ch = 1:2
        for color = 1:3
            PcolorGivenCh{s}{ch}{color} = sum(strcmpi(RDT.SubjectID, subjects{s}) & RDT.Channel == ch & RDT.Color == color)./...
                sum(strcmpi(RDT.SubjectID, subjects{s}) & RDT.Channel == ch & RDT.YesNo == 1);
        end

    end

Y(1,:) = horzcat(PcolorGivenCh{s}{2}{:});
Y(2,:) = horzcat(PcolorGivenCh{s}{1}{:});
Y = fliplr(Y);
Y(:,2:end) = fliplr(Y(:,2:end));

X = xpositions;
b{s} = bar(barAxes{s}, X, Y, 'stacked', 'LineWidth',1.25);
set(barAxes{s}, 'ColorOrder', [1 1 1; rgbMatrix_Red(end,:); rgbMatrix_Green(end,:)], 'linewidth', defaultLineWidth);
set(barAxes{s}, 'TickDir', 'out', 'box', 'off', 'TickLength', tickLength)
xlim(barAxes{s}, [0.1 0.9])
xticks(barAxes{s}, sort(xpositions(:)));

ylabel(barAxes{s}, 'Color response proportion', 'FontSize', 14)


if s == 2
xticklabels(barAxes{s}, {'543', '680'});
xlabel(barAxes{s}, 'Stimulus wavelength (nm)', 'FontSize', labelFontSize)
else
    xticklabels(barAxes{s}, {''})
end

end
for i = 1:numel(barAxes)
    barAxes{i}.FontSize = axFontSize;
    barAxes{i}.XColor = defaultColor; barAxes{i}.YColor = defaultColor;
    barAxes{i}.TickDir = tickDir;
end
    


%% P(expected response) vs heterogeneity

% get each subjects global ratio


inflectionProp= 0.5;
hetero = 1-(inflectionProp.^-1).*abs(propL-inflectionProp);
n = 2;
heteroEdges = [0 quantile(hetero, n) 1];
heteroEdges = round(2*heteroEdges,1)/2; % round to nearest 0.5


% for both subjects

seenIdx = RDT.YesNo == 1;
expectedIdx = RDT.Color == RDT.Channel;

[Nseen2, ~, ~, binHetero, binInt] = histcounts2(hetero(seenIdx), intensity(seenIdx), heteroEdges, intEdges);


uniqueBinInt = unique(binInt);
uniqueBinHetero = unique(binHetero);
uniqueBinInt(uniqueBinInt ==0) =  [];
uniqueBinHetero(uniqueBinHetero==0) = [];

heteroSeen = hetero(seenIdx);
intensitySeen = intensity(seenIdx);

for i = 1:numel(uniqueBinHetero)
    for j = 1:numel(uniqueBinInt)
        idx = binHetero == uniqueBinHetero(i) & binInt == uniqueBinInt(j);
        heteroMedian(i,j) = median(heteroSeen(idx));
        intensityMedian(i,j) = median(intensitySeen(idx));
    end
end


Nseen2 = Nseen2 .* (Nseen2 > minTrialsPerBin);
Nexpected2 = histcounts2(hetero(expectedIdx), intensity(expectedIdx), heteroEdges, intEdges);


Pexpected = Nexpected2./Nseen2;


Pexpected(~isfinite(Pexpected)) = nan;
intensityMidPoints = log2(intEdges(1:end-1)) + diff(log2(intEdges))/2;

figure; hold on;
markerEdgeColor = [0 0 0];
lineStyles = {':', '--', '-'};
markerFaceColor = [1 1 1; 0.5 0.5 0.5; 0 0 0];
maxMarkerSize = 12;
for i = 1:size(Pexpected,1) % for each heterogeneity

    color = i * (1/size(Pexpected,1));
    plot(intensityMedian(i,:), Pexpected(i,:), 'LineStyle', lineStyles{3} , 'Color', defaultColor, 'LineWidth',1.5);
    for j = 1:size(Pexpected,2)
        %markerSize = sqrt((Nseen2(i,j)./max(Nseen2(:)))*maxMarkerSize);
        markerSize = sqrt((Nseen2(i,j)./max(Nseen2(:))))*maxMarkerSize;

        if markerSize > 0
            plot(intensityMedian(i,j), Pexpected(i,j), 'Marker', 'o', 'LineWidth', 1.5,'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(i,:), 'MarkerSize', markerSize);
        else
        end
    end

end

set(gca, 'FontSize', axFontSize, 'xscale', 'log', 'yscale', 'log',  'XMinorTick', 'off', 'YMinorTick', 'off', 'Linewidth',2)
xlabel('Intensity (x threshold)')
ylabel('P(expected hue response)')
xlim([intEdges(2) intEdges(end)])
xticks(intEdges)

PexpectedTicks= 2.^(log2(0.0625):log2(2):log2(1));

yticks(PexpectedTicks);
YTickLab = cellfun(@(X) [num2str(X, '%0.3g') '%'], num2cell(100.*PexpectedTicks), 'UniformOutput', false);
yticklabels(YTickLab)


XTickLab = (cellfun(@(X) num2str(X, '%0.3g'), num2cell(round(intEdges,2)), 'UniformOutput', false));
XTickLab(1:2:end) = {''};
xticklabels(XTickLab);
ylim([0.0625/1.25 1.25])

axis square
set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5   9   9])
%set(gca, 'Position', [0.175 0.175 0.8 0.8])
set(gca, 'Position', [0.25 0.25 0.7 0.7])
set(gca, 'TickDir', 'out', 'TickLength', tickLength)
set(gcf, 'Color', 'w')
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters')
set(gca, 'FontSize', axFontSize)
xlabel('Intensity (x threshold)', 'FontSize', 14)
ylabel('P(expected color response) (%)', 'FontSize', 14)
set(gca, 'XTickLabelRotation', 0)

grid off;
%grid(gca);

%% Plot heterogeneity histogram
heteroSeen = hetero(seenIdx);

figure; hold on


for i = 1:numel(uniqueBinHetero)
hist{i} = histogram(heteroSeen(binHetero == uniqueBinHetero(i)), 'BinEdges', 0:0.05:1);
hist{i}.FaceAlpha = 1;
hist{i}.FaceColor = markerFaceColor(i,:);
end

set(gca, 'FontSize', axFontSize, 'XMinorTick', 'off', 'YMinorTick', 'off', 'Linewidth',2)
set(gca, 'Position', [0.25 0.25 0.7 0.7])
set(gcf, 'units', 'centimeters');
set(gcf, 'Position', [0.5   0.5   9   9]);
set(gca, 'TickDir', 'out', 'TickLength', tickLength)
set(gcf, 'Color', 'w')
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters')
set(gca, 'FontSize', axFontSize)
xlabel('Heterogeneity index', 'FontSize', 14);
ylabel('# Trials', 'FontSize', 14);
xlim([0 1]); ylim([0 450]);
set(gca, 'position', [0.25 0.25 0.7 0.3])