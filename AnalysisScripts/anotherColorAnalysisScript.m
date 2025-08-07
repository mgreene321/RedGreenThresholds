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
propThresh = nan(size(RDT,1),1);
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
propL = RDT.L./(RDT.L + RDT.M);
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
intensity = propThresh;
%intEdges = logspace(min(log10(propThresh(seen))), max(log10(propThresh(seen))), 4);
%intEdges = [0 1.5 30];
intEdges = [0 1.25 30];
%midPoints = intEdges(1:end-1) + diff(intEdges)/2;
%intensity = propThresh;
%propLEdges = [linspace(0, 0.95, 4) 1];
%propLEdges = [0 0.5 2/3 4/5 8/9 1];
propLEdges = linspace(0.25,1,6);
midPoints = propLEdges(1:end-1) + diff(propLEdges)/2;

[NpropL, ~,binsPropL] = histcounts(propL, 'BinEdges', propLEdges);

[Nintensity, ~, binsIntensity] = histcounts(intensity, 'BinEdges', intEdges);
 subjects = unique(ADT.SubjectID);
for s = 1:numel(subjects) % for each subject
    for ch = 1:numel(unique(RDT.Channel)) % for each channel
        for hue = 1:4 % for each hue response number

            trialIdx = strcmpi(RDT.SubjectID, subjects{s}) &...
                       RDT.Channel == ch & ...
                       RDT.Color == hue;

            trialsOfInterest(:,:,ch,hue,s) = trialIdx;

            % bin by propL and intensity
            N(:,:, ch, hue, s)= histcounts2(intensity(trialIdx), propL(trialIdx), intEdges, propLEdges);

        end
    end
end

% plot  680  responses for subj 2, eg
% resp680 = N(:,:,1,:,2);
% 
% resp680 = squeeze(resp680);
% resp680total = sum(resp680(:,:,1:3),3);
% resp680total(resp680total<5) = nan;
% Pred = resp680(:,:,1)./resp680total;
% Pgreen = resp680(:,:,2)./resp680total;
% Pachrom = resp680(:,:,3)./resp680total;
respMin = 5;
lineWidth = 2;
for s = 1:numel(subjects)
    figure;
    hold on
    plotCounter = 1;
    for ch = 1:2
        
        resp = N(:,:, ch,:,s);
        resp = squeeze(resp);
        respTotal = sum(resp(:,:,1:3),3);
        respTotal(respTotal<respMin) = nan;

        Pred = resp(:,:,1)./respTotal;
        Pgreen = resp(:,:,2)./respTotal;
        Pachrom = resp(:,:,3)./respTotal;

        for i = 1:numel(Nintensity)
            subplot(2, numel(Nintensity), plotCounter); hold on
       
            
            % get median propL for each propL bin

            for j = 1:numel(NpropL)
                medianPropL(i,j) = median(propL(trialsOfInterest(:,:,ch,hue,s) & binsIntensity == i & binsPropL == j));
            end


            xline(propLEdges, 'linestyle', '--', 'lineWidth',2, 'Color', [0 0 0])
            X = medianPropL(i,:);
            Y = [Pachrom(i,:)' Pgreen(i,:)' Pred(i,:)'];
            bar(X,Y, 'stacked');
            xlim([min(propLEdges) max(propLEdges)]);
            set(gca, 'ColorOrder', [1 1 1; 0 1 0; 1 0 0], 'linewidth', lineWidth);
            plotCounter = plotCounter+1;

        end


    end

end