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

subjectId = '20217R';
root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);
%RDT = RDT(strcmpi(RDT.SubjectID, subjectId),:);
% Load in aanalyzed table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
%ADT = ADT(strcmpi(ADT.SubjectID, subjectId),:);

% convert from intensity in au to "transmitted quanta" (quantal flux at
% photoreceptors)


% group table by folder
photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors = FillInPhotoreceptors(photoreceptors);
scanAngleDeg = 0.9;

%% Account for bleaching by projector background

whiteRGB = [82 90 128];

ND = 1;
projectorLum = (170.2 * 10^ND);

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
eyeParams.coneAperture = fspecial('gaussian', 11, gaussianSigma);

% Get cone quantum efficiencies at steady state bleached level, and
% projector retinal irradiance
[photoreceptorsBleached, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = ...
    projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);
wvl = SToWls(photoreceptorsBleached.nomogram.S);

luminanceContrast = nan(size(RDT,1), 1);

redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

% was SplineCMF
Red_energy = transpose(SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl));
Green_energy = transpose(SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl));

load T_xyz1931;
T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),photoreceptorsBleached.nomogram.S);
clear T_xyz1931 S_xyz1931



indices = getLocationIndicesInRawDataTable;
propThresh = nan(size(RDT,1),1);
for i = 1:size(RDT,1)
    temp = cellfun(@(x) find(x), indices, 'UniformOutput', false);
    ADT_idx = cellfun(@(x) ismember(i, x), temp);
    if RDT.Channel(i) == 1
        lambdaMax = 680;
        maxLaserPowerWatts = 0.4*unique(ADT.RedPowerWattsAt680nm(strcmpi(ADT.Folder, RDT.Folder(i))));
        energySpectrum = Red_energy;
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.RedThresholdPowerAU(ADT_idx);

    elseif RDT.Channel(i) == 2
        lambdaMax = 543;
        maxLaserPowerWatts = 0.4*unique(ADT.GreenPowerWattsAt543nm(strcmpi(ADT.Folder, RDT.Folder(i))));
        energySpectrum = Green_energy;
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.GreenThresholdPowerAU(ADT_idx);
    end

    %sigma = FWHM/(2*sqrt(2*log(2)));
    %energySpectrum = normpdf(wvl, lambdaMax, sigma);
    % 
    % valueInRealUnits = convertAOArbitraryToRealUnits(RDT.IntensityAU(i), maxLaserPowerWatts, lambdaMax, energySpectrum,photoreceptorsBleached.nomogram.S, scanAngleDeg);
    % 
    % 
    % radiance = valueInRealUnits.radianceWattsPerM2Srlambda;
    % luminance= 683.*dot(radiance, T_vLambda);
    % luminanceContrast(i,:) = 100.*luminance/projectorLum;%.* stimRenderTimeSec .* ISareaUm2;
end
intensity = propThresh;
%luminanceContrast = log10(luminanceContrast);
%incidentQuantaMax(incidentQuantaMax < 5) = 5;
%% Plot raw data
markerSize = 8;
fontSize = 20;
lineWidth = 3;
fontWeight = 'bold';
% alt way

% figure, hold on, set(gca, 'color', [1 1 1]./3); axis square
% colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0])
% set(gca, 'FontSize', fontSize, 'XColor', [0 0 0], 'YColor', [0 0 0], 'FontWeight', fontWeight, 'LineWidth', lineWidth, 'TickDir', 'out')
% xlim([0 1]);
% ylim([0 6]);
% xticks(0:(1/3):1);
% yticks(5:(2/3):9);
% xlabel('Proportion L');
% ylabel('Intensity (x threshold)')
% 
% R680 = RDT(RDT.Channel == 1 & strcmpi(RDT.SubjectID, subjectId),:);
% propL680 = R680.L./(R680.L + R680.M);
% intensity680 = intensity(RDT.Channel == 1);
% for i = 1:size(R680,1)
%     if R680.Color(i) == 1
%         scatter(propL680(i), intensity680(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha',1)
%     elseif R680.Color(i) == 2
%         scatter(propL680(i), intensity680(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 1 0], 'MarkerFaceAlpha', 1)
%     elseif R680.Color(i) == 3
%         scatter(propL680(i), intensity680(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1)
%     elseif R680.Color(i) == 4
%         %  scatter(propL680(i), transmittedQuanta680(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%     end
% 
% 
% end
% 
% figure, hold on, set(gca, 'color', [1 1 1]./3);axis square
% colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0])
% set(gca, 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k', 'FontWeight', fontWeight, 'LineWidth', lineWidth, 'TickDir', 'out')
% xlim([0 1]);
% ylim([0 6]);
% xticks(0:(1/3):1);
% yticks(5:(2/3):9);
% xlabel('Proportion L');
% ylabel('Intensity (x threshold)')
% 
% 
% R543 = RDT(RDT.Channel == 2 & strcmpi(RDT.SubjectID, subjectId),:);
% propL543 = R543.L./(R543.L + R543.M);
% intensity543 = intensity(RDT.Channel == 2);
% for i = 1:size(R543,1)
%     if R543.Color(i) == 1
%         scatter(propL543(i), intensity543(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha',1)
%     elseif R543.Color(i) == 2
%         scatter(propL543(i), intensity543(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 1 0], 'MarkerFaceAlpha',1)
%     elseif R543.Color(i) == 3
%         scatter(propL543(i), intensity543(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha',1)
%     elseif R543.Color(i) == 4
%         %scatter(propL543(i), transmittedQuanta543(i), markerSize^2, 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',0.5)
%     end
% 
% 
% end

propLEdges = linspace(0,1, 9);


minIntensity= min(intensity(RDT.YesNo==1));
maxIntensity= max(intensity(RDT.YesNo==1));

intEdges = 10.^(log10(minIntensity):log10(2):log10(maxIntensity));


%% colormaps
numColorMapRows = 256;

rgbMatrix_Red = [linspace(0,1, 256)', linspace(0, 0.25, 256)', linspace(0,0.25, 256)'];
rgbMatrix_Green = [linspace(0, 0.33, 256)', linspace(0,0.9, 256)', linspace(0, 0.33, 256)'];
rgbMatrix_Yellow = [linspace(0.5, 1, 256)' linspace(0.5, 1, 256)' 0.5*(ones(256,1))];

propL = RDT.L./(RDT.L + RDT.M);
for ch = 1:2
    for col = 1:4
        idx = strcmpi(RDT.SubjectID, subjectId) & RDT.Channel == ch & RDT.Color == col;
        N{ch, col} = histcounts2(propL(idx), intensity(idx),propLEdges, intEdges);
    end
end

minTrialsPerBin =10;
for ch = 1:2

    % Distribution of all seen trials
    %     figure; hold on;
    %     imagesc(rot90(sum(cat(3, N{ch,1:3}),3)));
    Nseen{ch} = sum(cat(3, N{ch, 1:3}),3);
    Nseen{ch} = Nseen{ch} .* (Nseen{ch} >= minTrialsPerBin);
    NexpectedColor{ch} = N{ch,ch};
    temp = NexpectedColor{ch}./Nseen{ch};
    % handle nan and inf

    temp = log2(temp);
    temp(temp == -Inf) = -10;
    temp(isnan(temp) | isinf(temp)) = nan;
    temp = rot90(temp);
    temp = flipud(temp);
   

    % figure; hold on; set(gca, 'ydir', 'reverse'); axis square;
    % 
    % set(gca, 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k', 'FontWeight', fontWeight, 'LineWidth', lineWidth)
    % imagesc(rot90(temp));
    % axis square

    %
       %
    numRows = size(temp,1);
    numCols = size(temp,2);
    xVals = 0:1:numCols;
    yVals = 0:1:numRows;

[xGrid, yGrid] = meshgrid(xVals, yVals);

figure; set(gcf, 'renderer', 'painters', 'paperpositionmode', 'auto', 'color', 'w');


  linesPerSquareEdge = 4;
numBoxes = max([numCols numRows]); % Make it square for the largest dimension of the data matrix
        plot([zeros(1,linesPerSquareEdge*(numBoxes+1)); linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))],...
            [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))], 'k-'); hold on;
        plot([zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes; linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1))],...
            [linspace(0,numBoxes,linesPerSquareEdge*(numBoxes+1)); zeros(1,linesPerSquareEdge*(numBoxes+1))+numBoxes], 'k-'); hold on;

        for r = 1:numRows
            for c = 1:numCols
                if isfinite(temp(r,c))
                    patch([xGrid(r,c:c+1) fliplr(xGrid(r+1,c:c+1))], [yGrid(r,c:c+1) fliplr(yGrid(r+1,c:c+1))], temp(r,c), 'EdgeColor', 'w', 'linewidth', 1.25);
                end
            end
        end




    %
    %


     xlim([0 numel(propLEdges)-1]);
     ylim([0 numel(intEdges)-1]);
    % xx = 0:1:(numel(propLEdges)-1);
    % yy = 1:1:(numel(intEdges));
    % yticks(yy(1:2:end));
    % xticks(xx(1:2:end));

    propLTicks = 0+mod(numel(propLEdges)-1, 2):2: numel(propLEdges) -mod(numel(propLEdges)-1, 2);
    intTicks = 0+mod(numel(intEdges)-1, 2):2: numel(intEdges) -mod(numel(intEdges)-1, 2);

    xticks(propLTicks);
    yticks(intTicks);


    xTickLabs =  cellfun(@(X) num2str(X, '%0.3g'), num2cell(propLEdges(propLTicks+1)), 'UniformOutput', false);
    yTickLabs =  cellfun(@(X) num2str(X, '%0.3g'), num2cell(round(intEdges(intTicks+1),2)), 'UniformOutput', false);

    xticklabels(xTickLabs);
    yticklabels(yTickLabs);

    % yticklabels(cellfun(@(X) num2str(X, '%0.2f'),num2cell(intEdges(1:2:end)), 'UniformOutput', false))
    % %yticklabels(fliplr({'5', '5.67', '6.33', '7', '7.67', '8.33', '9'}))
    % xticklabels(fliplr(cellfun(@(X) num2str(X, '%0.2f'), num2cell(propLEdges(1:2:end)), 'UniformOutput', false)))
    %xticklabels({'0', '0.33', '0.67', '1'})
    xlabel('Proportion L');
    ylabel('Intensity (x threshold)');


    if ch == 1

        c = rgbMatrix_Red;

       % c(end+1,:) =  [1 1 1];
        colormap(gca, c);
       
        cb = colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0], 'Location', 'northoutside');
        %caxis([0 1.01]);

    elseif ch == 2
        c = rgbMatrix_Green;
        %c(end+1,:) = [1 1 1];
        colormap(gca, c);

        cb = colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0], 'Location', 'northoutside');
        %caxis([0 1.01]);
    end
    % 
    % xline(0.5:1:(numel(propLEdges)-0.5), 'k', 'LineWidth', 2)
    % yline(0.5:1:(numel(intEdges)-0.5), 'k', 'LineWidth', 2)

    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [0.5+ch/2   0.5    7    7])

    set(gca, 'position', [0.2 0.2 0.6 0.6]);
    set(gca, 'FontSize', fontSize)
    %set(gca, 'ydir', 'reverse')
    box on;
   PexpectedTicks= 2.^(log2(0.0625):log2(2):log2(1));
    clim([log2(0.0625/1.25) log2(1.25)]);
    cb.Ticks = [-4:1:0];
    cb.TickLabels = cellfun(@(X) [num2str(X, '%0.3g') '%'], num2cell(100.*PexpectedTicks), 'UniformOutput', false);
    set(gca, 'LineWidth', 2)


end

% P(expected response) vs heterogeneity

% get each subjects global ratio
subjects = unique(RDT.SubjectID);

%inflectionProp = 0.5;
% for s = 1:numel(subjects)
%     [nL(s), nM(s), nS(s)] = getNumberLMSCones(subjects{s});
%     pL(s) = nL(s)./(nL(s) + nM(s));
% 
%     inflectionProp(strcmpi(RDT.SubjectID, subjects{s}),:) = pL(s);
% end


inflectionProp= 0.5;
hetero = 1-(inflectionProp.^-1).*abs(propL-inflectionProp);
n = 2;
heteroEdges = [0 quantile(hetero, n) 1];
heteroEdges = round(2*heteroEdges,1)/2; % round to nearest 0.5

seenIdx = RDT.YesNo == 1 & strcmpi(RDT.SubjectID, subjectId);
expectedIdx = RDT.Color == RDT.Channel & strcmpi(RDT.SubjectID, subjectId);

[Nseen2, ~, ~, binInt, binHetero] = histcounts2(hetero(seenIdx), intensity(seenIdx), heteroEdges, intEdges);
Nexpected2 = histcounts2(hetero(expectedIdx), intensity(expectedIdx), heteroEdges, intEdges);
Nseen2 = Nseen2 .* (Nseen2 > minTrialsPerBin);

% get median values in each bin

Pexpected = Nexpected2./Nseen2;
Pexpected(isnan(Pexpected) | isinf(Pexpected)) = nan;
% Pexpected(isnan(Pexpected) | isinf(Pexpected)) = 1.01;


figure; hold on; set(gca, 'ydir', 'reverse'); axis square;

    set(gca, 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k', 'FontWeight', fontWeight, 'LineWidth', lineWidth)
    imagesc(rot90(Pexpected));

 

    axis square
xlim([0.5 numel(heteroEdges)-0.5]);
    ylim([0.5 numel(intEdges)-0.5]);
    xx = 0.5:1:(numel(heteroEdges)-0.5);
    yy = 0.5:1:(numel(intEdges)-0.5);
    yticks(yy(1:2:end));
    xticks(xx(1:2:end));
    yticklabels(cellfun(@(X) num2str(X, '%0.3g'), fliplr(num2cell(intEdges(2:2:end))), 'UniformOutput', false))
    %yticklabels(fliplr({'5', '5.67', '6.33', '7', '7.67', '8.33', '9'}))
    xticklabels(fliplr(cellfun(@(X) num2str(X, '%0.3g'), fliplr(num2cell(heteroEdges(1:2:end))), 'UniformOutput', false)))
    %xticklabels({'0', '0.33', '0.67', '1'})
    xlabel('Heterogeneity');
    ylabel('Intensity (x threshold)');
    c = rgbMatrix_Yellow;
        c(end+1,:) =  [1 1 1];
        colormap(gca, c);
        colorbar;
        colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0])
        caxis([0 1.01]);

        set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [0.5+ch/2   0.5    7    7])

    set(gca, 'position', [0.2 0.2 0.6 0.6]);
%  yline(0.5:1:(numel(intEdges)-0.5), 'k', 'LineWidth', 2)
% xline(0.5:1:(numel(heteroEdges)-0.5), 'k', 'LineWidth', 2)


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
Pexpected(isnan(Pexpected) | isinf(Pexpected)) = 1.01;


figure; hold on; set(gca, 'ydir', 'reverse'); axis square;

    set(gca, 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k', 'FontWeight', fontWeight, 'LineWidth', lineWidth)
    imagesc(rot90(Pexpected));

    axis square
xlim([0.5 numel(heteroEdges)-0.5]);
    ylim([0.5 numel(intEdges)-0.5]);
    xx = 0.5:1:(numel(heteroEdges)-0.5);
    yy = 0.5:1:(numel(intEdges)-0.5);
    yticks(yy(1:2:end));
    xticks(xx(1:2:end));
    yticklabels(cellfun(@(X) num2str(X, '%0.3g'), fliplr(num2cell(intEdges(2:2:end))), 'UniformOutput', false))
    %yticklabels(fliplr({'5', '5.67', '6.33', '7', '7.67', '8.33', '9'}))
    xticklabels(fliplr(cellfun(@(X) num2str(X, '%0.3g'), fliplr(num2cell(heteroEdges(1:2:end))), 'UniformOutput', false)))
    %xticklabels({'0', '0.33', '0.67', '1'})
    xlabel('Heterogeneity');
    ylabel('Intensity (x threshold)');
    c = rgbMatrix_Yellow;
        c(end+1,:) =  [1 1 1];
        colormap(gca, c);
        colorbar;
        colorbar('LineWidth', lineWidth, 'FontSize', fontSize, 'FontWeight', fontWeight, 'TickDirection', 'out', 'Color', [0 0 0])
        caxis([0 1.01]);

        set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [0.5+ch/2   0.5    7    7])

    set(gca, 'position', [0.2 0.2 0.6 0.6]);


%% 

Pexpected2 = Pexpected;
Pexpected2(Pexpected2 == 1.01) = nan;
intensityMidPoints = log2(intEdges(1:end-1)) + diff(log2(intEdges))/2;

figure; hold on;
markerEdgeColor = [0 0 0];
lineStyles = {':', '--', '-'};
markerFaceColor = [1 1 1; 0.5 0.5 0.5; 0 0 0];
maxMarkerSize = 600;
for i = 1:size(Pexpected2,1) % for each heterogeneity

    color = i * (1/size(Pexpected2,1));
    plot(intensityMedian(i,:), Pexpected2(i,:), 'LineStyle', lineStyles{3} , 'Color', 'k', 'LineWidth',2);
for j = 1:size(Pexpected,2)
    markerSize = sqrt((Nseen2(i,j)./max(Nseen2(:)))*maxMarkerSize);

    if markerSize > 0
    plot(intensityMedian(i,j), Pexpected2(i,j), 'Marker', 'o', 'LineWidth', 1.5,'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor', markerFaceColor(i,:), 'MarkerSize', markerSize);
    else
    end
end

end

set(gca, 'FontSize', fontSize, 'xscale', 'log', 'yscale', 'log',  'XMinorTick', 'off', 'YMinorTick', 'off', 'Linewidth',2)
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
set(gcf, 'Units', 'inches')
set(gcf, 'Position', [0.5   0.5    7    7])

set(gca, 'Position', [0.175 0.175 0.8 0.8])
grid off;
grid(gca);