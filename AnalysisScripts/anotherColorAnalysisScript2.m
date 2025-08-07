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
subjects = unique(ADT.SubjectID);
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

intensity = propThresh;
indices = getLocationIndicesInRawDataTable;

propL = RDT.L./(RDT.L + RDT.M);

%globalPropL = 

inflectionProp = 0.5;
hetero = 1-(inflectionProp^-1).*abs(propL-inflectionProp);


heteroEdges = linspace(0,1,4);
%intEdges = logspace(log10(0.4), log10(30), 4);
intEdges = [min(intensity) median(intensity) max(intensity)];
heteroMidPoints = heteroEdges(1:end-1) + diff(heteroEdges)/2;
% bin by het idx & intensity
intensity = propThresh;

% lum contrast
%intensity = luminanceContrast;
%intEdges = logspace(1,2.8,10);
%%
%[NN, ~, ~, binsHetero, binsIntensity] = histcounts2(hetero, intensity, heteroEdges, intEdges);
[~,~, binsHetero] = histcounts(hetero, heteroEdges);
[~, ~, binsIntensity] = histcounts(intensity, intEdges);

uniqueBinsHetero = unique(binsHetero);
uniqueBinsIntensity = unique(binsIntensity);

% for each subject
for s = 1:numel(subjects)
    % for each channel
    for ch = 1:2
        for color = 1:3
            % for each hetero bin
            for h = 1:numel(uniqueBinsHetero)
                % for each intensity bin
                for i = 1:numel(uniqueBinsIntensity)
                    idx = strcmpi(RDT.SubjectID, subjects{s}) &...
                        RDT.Channel == ch &...
                        RDT.Color == color &...
                        binsHetero == uniqueBinsHetero(h) &...
                        binsIntensity == uniqueBinsIntensity(i);

                    Presp{s}{ch}{color}(h,i) = sum(idx)./sum(strcmpi(RDT.SubjectID, subjects{s}) &...
                        RDT.Channel == ch &...
                        binsHetero == uniqueBinsHetero(h) &...
                        binsIntensity == uniqueBinsIntensity(i));

                    if sum(strcmpi(RDT.SubjectID, subjects{s}) &...
                        RDT.Channel == ch &...
                        binsHetero == uniqueBinsHetero(h) &...
                        binsIntensity == uniqueBinsIntensity(i)) < 10
                        Presp{s}{ch}{color}(h,i) = nan;
                    end

                end

            end

        end

    end
end

% prob of veridical response

for s = 1:numel(subjects)
    for h = 1:numel(uniqueBinsHetero)
        for i = 1:numel(uniqueBinsIntensity)

            idx = strcmpi(RDT.SubjectID, subjects{s}) &...
                binsHetero == uniqueBinsHetero(h) &...
                binsIntensity == uniqueBinsIntensity(i) &...
                RDT.Channel == RDT.Color;


            Pver{s}(h,i) = sum(idx)./sum(strcmpi(RDT.SubjectID, subjects{s}) &...
                binsHetero == uniqueBinsHetero(h) &...
                binsIntensity == uniqueBinsIntensity(i));

            if sum(strcmpi(RDT.SubjectID, subjects{s}) &...
                    binsHetero == uniqueBinsHetero(h) &...
                    binsIntensity == uniqueBinsIntensity(i)) < 10

                Pver{s}(h,i) = nan;
            end

        end
    end

end

% plot for each loc

[~,q] = iqr(propThresh(RDT.YesNo==1));

maxInt = 32;
minInt = 0;
minTrials = 10;
for loc = 1:size(ADT,1)
    idx =   indices{loc} &...
       propThresh > minInt & propThresh <maxInt &...
        RDT.Channel == RDT.Color;

    Pver(loc) = sum(idx)./sum(indices{loc} &...
         propThresh > minInt & propThresh <maxInt & RDT.YesNo == 1);



    if sum(indices{loc} &...
          propThresh > minInt & propThresh <maxInt & RDT.YesNo == 1) < minTrials
        Pver(loc) = nan;
    end


end


for loc = 1:size(ADT,1)
    idx =   indices{loc} &...
        propThresh > minInt & propThresh < maxInt &...
        RDT.Channel == 1 & RDT.Color== 1;

    PverRed(loc) = sum(idx)./sum(indices{loc} &...
         propThresh > minInt & propThresh < maxInt & RDT.Channel == 1 & RDT.YesNo == 1);


    
    if sum(indices{loc} &...
             propThresh > minInt & propThresh < maxInt & RDT.Channel == 1 & RDT.YesNo == 1) < minTrials
        PverRed(loc) = nan;
    end


end


for loc = 1:size(ADT,1)
    idx =   indices{loc} &...
        propThresh <1.25 &...
        RDT.Channel == 2 & RDT.Color== 2;

    PverGreen(loc) = sum(idx)./sum(indices{loc} &...
         propThresh > minInt & propThresh < maxInt  & RDT.Channel == 2 & RDT.YesNo == 1);
    
    if sum(indices{loc} &...
             propThresh > minInt & propThresh < maxInt & RDT.Channel == 2 & RDT.YesNo == 1) < minTrials
        PverGreen(loc) = nan;
    end


end
