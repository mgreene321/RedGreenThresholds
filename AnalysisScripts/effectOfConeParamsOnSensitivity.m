pL = 0.9;

%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
pupilDiamMm = 6.5;
pupilDiamCm = pupilDiamMm*(10^-1);
eyeLengthMm = 16.7;
ISdiameterUm = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)
stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
stimSideLengthPix = 21;
pixPerDeg = 560;
stimSideLengthDeg = stimSideLengthPix/pixPerDeg;
stimSideLengthUm = DegreesToRetinalMM(eyeLengthMm, stimSideLengthDeg) * 1e3;
stimAreaUm2 = stimSideLengthUm^2;
ISareaUm2 = pi.*(ISdiameterUm/2).^2;

ppd =560;
coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneDiamDeg = coneDiamPixels./ppd;
coneDiamUm = DegreesToRetinalMM(coneDiamDeg, eyeLengthMm)*1e3;
coneAreaUm2 = pi*(coneDiamUm/2)^2;
%% Utility calculations
scanAngleDeg = 0.9;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pi*((pupilDiamCm/2)^2);
eyeLengthCm = eyeLengthMm*(10^-1);
%%

eccDeg = linspace(0.1, 5, 100); % deg
stimSize = 2.25/60; % deg

photoreceptors_default = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors_default.macularPigmentDensity.source = 'CIE';
photoreceptors_default.lensDensity.source = 'CIE';

%bleaching

%% Account for bleaching by projector background

% Make cone aperture
coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 11;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma);

whiteRGB = [82 90 128];
%ND = -log10(20.39/170.2);
ND = 1;
projectorLum = (170.2 * 10^ND);
%projectorLum = 100;

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
%eyeParams.ISdiameterUm = coneDiamUm;
eyeParams.coneAperture = coneAperture;

ppd = 560; % pixels per degree
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;

%photoreceptorsBleached = projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);

[fractionBleachedFromIsom, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = ...
    projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);


photoreceptors_default.fractionPigmentBleached.value = fractionBleachedFromIsom;


photoreceptors = FillInPhotoreceptors(photoreceptors_default);


wls = SToWls(photoreceptors.nomogram.S);
% alter macular pigment density
mpodFig = figure; hold on
title('Vary MPOD')
xlabel('Eccentricity (deg)')
ylabel('S_{680}/S_{543}')
mpodAx = gca;

vpodFig = figure; hold on
title('Vary VPOD (axial density)')
xlabel('Eccentricity (deg)')
ylabel('S_{680}/S_{543}')
vpodAx = gca;
for ecc = 1:numel(eccDeg)

    % macular pigment
    mpod = computePeakMacularPigmentDensity(eccDeg(ecc), stimSize);

    ph_mpod{ecc} = photoreceptors;
    ph_mpod{ecc}.macularPigmentDensity.density = mpod .* ...
        ph_mpod{ecc}.macularPigmentDensity.density./max(ph_mpod{ecc}.macularPigmentDensity.density);

    ph_mpod{ecc}.macularPigmentDensity.transmittance = 10.^-ph_mpod{ecc}.macularPigmentDensity.density;

    ph_mpod{ecc}.preReceptoral.transmittance =  ph_mpod{ecc}.lensDensity.transmittance .* ...
        ph_mpod{ecc}.macularPigmentDensity.transmittance;

    ph_mpod{ecc}.effectiveAbsorptance =  ph_mpod{ecc}.absorptance .* ...
        (ones(size(ph_mpod{ecc}.absorptance,1),1)*ph_mpod{ecc}.preReceptoral.transmittance);

    for i = 1:size(ph_mpod{ecc}.effectiveAbsorptance,1)
        ph_mpod{ecc}.isomerizationAbsorptance(i,:) =  ph_mpod{ecc}.quantalEfficiency.value(i) * ...
            ph_mpod{ecc}.effectiveAbsorptance(i,:);
    end

    ph_mpod{ecc}.energyFundamentals = EnergyToQuanta(S, ph_mpod{ecc}.isomerizationAbsorptance')';
    mx = max(ph_mpod{ecc}.energyFundamentals,[],2);
    ph_mpod{ecc}.energyFundamentals = diag(1./mx)*ph_mpod{ecc}.energyFundamentals;

    %% Compute normalized quantal sensitivities (aka cone fundamentals in quantal units)
    ph_mpod{ecc}.quantalFundamentals = ph_mpod{ecc}.isomerizationAbsorptance;
    mx = max(ph_mpod{ecc}.quantalFundamentals,[],2);
    ph_mpod{ecc}.quantalFundamentals = diag(1./mx)*ph_mpod{ecc}.quantalFundamentals;

    S680_mpod = pL.*ph_mpod{ecc}.isomerizationAbsorptance(1,wls == 680) + (1-pL).*ph_mpod{ecc}.isomerizationAbsorptance(2,wls == 680);
    S543_mpod =  pL.*ph_mpod{ecc}.isomerizationAbsorptance(1,wls == 543) + (1-pL).*ph_mpod{ecc}.isomerizationAbsorptance(2,wls == 543);


    plot(mpodAx, eccDeg(ecc), S680_mpod/S543_mpod, 'ko');

    %%%%
    ph_vpod{ecc} = photoreceptors_default;
    ph_vpod{ecc}.specificDensity.source = 'None';
    ph_vpod{ecc}.OSlength.source = 'None';
    vpod = computePeakVisualPigmentDensity(eccDeg(ecc), stimSize);

    ph_vpod{ecc}.axialDensity.value = vpod;
    ph_vpod{ecc} = FillInPhotoreceptors(ph_vpod{ecc});

    S680_vpod = pL.*ph_vpod{ecc}.isomerizationAbsorptance(1,wls == 680) + (1-pL).*ph_vpod{ecc}.isomerizationAbsorptance(2,wls == 680);
    S543_vpod =  pL.*ph_vpod{ecc}.isomerizationAbsorptance(1,wls == 543) + (1-pL).*ph_vpod{ecc}.isomerizationAbsorptance(2,wls == 543);
    plot(vpodAx, eccDeg(ecc), S680_vpod/S543_vpod, 'ko');


end

% now try just altering OSlength

figure; hold on
title('Vary OS length')
xlabel('OS length (\mum)')
ylabel('S_{680}/S_{543}')
osAx = gca;
OSlengths = (10:30).*[1 1 1]';

for o = 1:size(OSlengths,2)
    ph_os{o} = photoreceptors_default;
    ph_os{o}.OSlength.value = OSlengths(:,o);

    ph_os{o} = FillInPhotoreceptors(ph_os{o});
    S680_os = pL.*ph_os{o}.isomerizationAbsorptance(1,wls == 680) + (1-pL).*ph_os{o}.isomerizationAbsorptance(2,wls == 680);
    S543_os =  pL.*ph_os{o}.isomerizationAbsorptance(1,wls == 543) + (1-pL).*ph_os{o}.isomerizationAbsorptance(2,wls == 543);
    plot(osAx, unique(OSlengths(:,o)), S680_os/S543_os, 'ro');

end

% % alter just ISdiameter
% 
% figure; hold on
% title('Vary IS diameter')
% xlabel('OS length (\mum)')
% ylabel('S_{680}/S_{543}')
% isAx = gca;
% ISdiams = (linspace(2,5,50)).*[1 1 1]';
% 
% for i = 1:size(ISdiams,2)
%     ph_is{i} = photoreceptors_default;
%     ph_is{i}.ISdiameter.value = ISdiams(:,i);
% 
%     ph_is{i} = FillInPhotoreceptors(ph_is{i});
%     S680_is = ph_is{i}.isomerizationAbsorptance(1,wls == 680) + ph_is{i}.isomerizationAbsorptance(2,wls == 680);
%     S543_is =  ph_is{i}.isomerizationAbsorptance(1,wls == 543) + ph_is{i}.isomerizationAbsorptance(2,wls == 543);
%     plot(isAx, unique(ISdiams(:,i)), S680_is/S543_is, 'ro');
% 
% end