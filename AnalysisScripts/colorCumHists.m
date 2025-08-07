% Load analyzed data table
scanAngleDeg = 0.9;
 S = [380 1 401];
 wvl = SToWls(S);
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
%Load analyzed data table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);

% load raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);

% compute luminance contrast on each trial

luminanceContrast = nan(size(RDT,1), 1);

redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

% was SplineCMF
Red_energy = transpose(SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl));
Green_energy = transpose(SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl));

load T_xyz1931;
T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
clear T_xyz1931 S_xyz1931

ND = 1;
projectorLum = (170.2 * 10^ND);

indices = getLocationIndicesInRawDataTable;
for i = 1:size(RDT,1)

    % find correspodning data point in ADT
    temp = cellfun(@(x) find(x), indices, 'UniformOutput', false);
    ADT_idx = cellfun(@(x) ismember(i, x), temp);
    
    if RDT.Channel(i) == 1
        lambdaMax = 680;
        maxLaserPowerWatts = 0.4*unique(ADT.RedPowerWattsAt680nm(ADT_idx));    
        energySpectrum = Red_energy;
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.RedThresholdPowerAU(ADT_idx);
        
    elseif RDT.Channel(i) == 2
        lambdaMax = 543;
        maxLaserPowerWatts = 0.4*unique(ADT.GreenPowerWattsAt543nm(ADT_idx));
        energySpectrum = Green_energy;
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.GreenThresholdPowerAU(ADT_idx);
    end
    
    %sigma = FWHM/(2*sqrt(2*log(2)));
    %energySpectrum = normpdf(wvl, lambdaMax, sigma);
    
    valueInRealUnits = convertAOArbitraryToRealUnits(RDT.IntensityAU(i), maxLaserPowerWatts, lambdaMax, energySpectrum,S, scanAngleDeg);
   
    
    radiance = valueInRealUnits.radianceWattsPerM2Srlambda;
    luminance= 683.*dot(radiance, T_vLambda);
    luminanceContrast(i,:) = 100.*luminance/projectorLum;%.* stimRenderTimeSec .* ISareaUm2;

end


indices = getLocationIndicesInRawDataTable;
propL = ADT.L./(ADT.L + ADT.M);
% 
% normPropL = nan(size(ADT,1),1);
% normPropL(strcmpi(ADT.SubjectID, '10001R')) = globalPropLWill;
% normPropL(strcmpi(ADT.SubjectID, '20217R')) = globalPropLMax;
% 
% normPropL = propL./normPropL;

% group by propL
%propLEdges =  [0 0.1 0.4 0.6 0.9 1];
%propLEdges = 0:(1/3):1;
%propLEdges = [0 1/3 2/3 (2/3 + 1/6 + 1/12) 3/3];

propLEdges = [0 0.3 0.6 0.9 1];
[NpropL, ~,binsPropL] = histcounts(propL, 'BinEdges', propLEdges);

% for luminance contrast
linearEdges = 3:10;
%intEdges = 2.^linearEdges;
linearMidPoints = linearEdges(1:end-1) + diff(linearEdges)/2;
midPoints = 2.^linearMidPoints;
%intensity = luminanceContrast;
intensity = propThresh;

% for threshold
intEdges = logspace(min(log10(propThresh(seen))), max(log10(propThresh(seen))), 7);
midPoints = intEdges(1:end-1) + diff(intEdges)/2;
%intensity = propThresh;

[Nintensity, ~, binsIntensity] = histcounts(intensity, 'BinEdges', intEdges);


% go through and get luminanceContrast of trials in each bin

% color response indices

red680 = RDT.Channel == 1 & RDT.Color == 1;
green680 = RDT.Channel == 1 & RDT.Color == 2;
achrom680 = RDT.Channel == 1 & RDT.Color == 3;
seen680 = RDT.Channel == 1 & RDT.YesNo ==1;

red543 = RDT.Channel == 2 & RDT.Color == 1;
green543 = RDT.Channel == 2 & RDT.Color == 2;
achrom543 = RDT.Channel == 2 & RDT.Color == 3;
seen543 = RDT.Channel == 2 & RDT.YesNo == 1;

seen = seen543 | seen680;
% get subject indices

subjects = unique(ADT.SubjectID);
subjIdx = cell(numel(subjects),1);
for s = 1:numel(subjects)
    subjIdx{s} = strcmpi(RDT.SubjectID, subjects{s});
end




[binSeen680, binRed680, binGreen680, binAchrom680,...
    binSeen543, binRed543, binGreen543, binAchrom543,...
    Nseen680, Nred680, Ngreen680, Nachrom680,...
    Nseen543, Nred543, Ngreen543, Nachrom543] = deal(cell(max(binsPropL), numel(subjects)));


for bin = 1:numel(NpropL)
    ADTnums = find(binsPropL==bin);
    temp = indices(ADTnums); % cell array containing RDT indices for each element in ADT indicated by ADTnums
    temp2 = cell2mat(temp);
    temp2 = reshape(temp2, [], size(temp,1));
    RDTindices = any(temp2, 2);
    if isempty(RDTindices)
        RDTindices = zeros(size(RDT,1),1);
    else
    end
    
   for s  = 1:numel(subjects)
   %binContrasts{bin,:} = luminanceContrast(RDTindices,:);
   binSeen680{bin,s} = intensity(RDT   indices & seen680 & subjIdx{s});
   binRed680{bin,s} = intensity(RDTindices & red680 & subjIdx{s});
   binGreen680{bin,s} = intensity(RDTindices & green680 & subjIdx{s});
   binAchrom680{bin,s} =intensity(RDTindices & achrom680 & subjIdx{s});

   binSeen543{bin,s} = intensity(RDTindices & seen543 & subjIdx{s});
   binRed543{bin,s} =intensity(RDTindices & red543 & subjIdx{s});
   binGreen543{bin,s} = intensity(RDTindices & green543 & subjIdx{s});
   binAchrom543{bin,s} = intensity(RDTindices & achrom543 & subjIdx{s});

   Nseen680{bin,s} = histcounts(binSeen680{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Nred680{bin,s} = histcounts(binRed680{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Ngreen680{bin,s} = histcounts(binGreen680{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Nachrom680{bin,s} = histcounts(binAchrom680{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');

   Nseen543{bin,s} = histcounts(binSeen543{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Nred543{bin,s} = histcounts(binRed543{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Ngreen543{bin,s} = histcounts(binGreen543{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');
   Nachrom543{bin,s} = histcounts(binAchrom543{bin,s}, 'BinEdges', intEdges, 'normalization', 'cumcount');

   Nseen680{bin,s}(Nseen680{bin,s} <10) = nan;
   Nseen543{bin,s}(Nseen543{bin,s} < 10) = nan;

 

   end
end

%plot
lineStyles = {'-', '--'};
lineWidth = 2;
markerSize = 8;
fontSize = 14;
for bin = 1:numel(NpropL)
    figure; hold on; grid on; 
     xlabel('Intensity (x threshold)');
     ylabel('P(color response | intensity <= x)');
    %xline(intEdges, 'k', 'LineWidth', lineWidth);
    for s  = 1:numel(subjects)
        plot(intEdges(2:end), Nred680{bin,s}./Nseen680{bin,s}, 'ro',  'MarkerFaceColor', 'r','LineStyle',  lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
        plot(intEdges(2:end), Ngreen680{bin,s}./Nseen680{bin,s}, 'go',  'MarkerfaceColor', 'g', 'LineStyle', lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
        plot(intEdges(2:end), Nachrom680{bin,s}./Nseen680{bin,s}, 'MarkerFaceColor', (2/3).*[1 1 1], 'Color', (2/3).*[1 1 1],'Marker', 'o', 'MarkerEdgeColor', (2/3).*[1 1 1],  'LineStyle', lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
    end
    xlim([min(intEdges) max(intEdges)]);
   ylim([0 1]);
   set(gca, 'xscale', 'log', 'LineWidth', lineWidth, 'xminorgrid', 'off', 'yminorgrid', 'off', 'fontsize', fontSize);
   xticks(intEdges);
   yticks(0:0.25:1);
   xticklabels(cellfun(@(x) num2str(x, '%0.2f'), num2cell(intEdges), 'UniformOutput', false));
  
end

for bin = 1:numel(NpropL)
    figure; hold on; grid on
    xlabel('Intensity (x threshold)');
    ylabel('P(color response | intensity <= x)');
       %xline(intEdges, 'k', 'LineWidth', lineWidth);
    for s  = 1:numel(subjects)
        plot(intEdges(2:end), Nred543{bin,s}./Nseen543{bin,s}, 'ro', 'MarkerFaceColor', 'r', 'LineStyle',  lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
        plot(intEdges(2:end), Ngreen543{bin,s}./Nseen543{bin,s}, 'go', 'MarkerFaceColor', 'g', 'LineStyle', lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
        plot(intEdges(2:end), Nachrom543{bin,s}./Nseen543{bin,s}, 'MarkerFaceColor',(2/3).*[1 1 1], 'Color', (2/3).*[1 1 1], 'Marker', 'o', 'MarkerEdgeColor', (2/3).*[1 1 1], 'LineStyle', lineStyles{s}, 'lineWidth', lineWidth, 'MarkerSize', markerSize);
    end
     xlim([min(intEdges) max(intEdges)]);
    ylim([0 1]);
      set(gca, 'xscale', 'log', 'LineWidth', lineWidth, 'xminorgrid', 'off', 'yminorgrid', 'off', 'fontsize', fontSize);
   xticks(intEdges);
     yticks(0:0.25:1);
     xticklabels(cellfun(@(x) num2str(x, '%0.2f'), num2cell(intEdges), 'UniformOutput', false));

end

