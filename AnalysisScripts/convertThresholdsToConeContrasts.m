function coneContrasts = convertThresholdsToConeContrasts

S = [380 1 401];
wvl = SToWls(S);
%% Physical constants

h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
%% Eye constants

pupilDiamMm = 6.5;
eyeLengthMm = 16.7;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);
photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Resolution of imaging system
ppd = 560; % pixels per degree
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;

%% Create cone aperture function, whose volume = 1

coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 712;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma); % Gaussian cone aperture
coneAperture = coneAperture./max(coneAperture(:));
whiteRGB = [82 90 128];
%ND = -log10(20.39/170.2);
ND = 1;
projectorLum = 170.2 * 10^ND;

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
eyeParams.coneAperture = fspecial('gaussian', 11, gaussianSigma); 

[photoreceptorsBleached, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);
% Background retinal power per pixel
bkgScalar = sum(projectorRetIrradianceQuantaPerLambdaUm2Sec(:));
normBkgRetIrrad = projectorRetIrradianceQuantaPerLambdaUm2Sec./bkgScalar;


L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
M = photoreceptorsBleached.isomerizationAbsorptance(2,:);
S = photoreceptorsBleached.isomerizationAbsorptance(3,:);

%% Create cone aperture function, whose volume = 1

coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 712;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma); % Gaussian cone aperture
coneAperture = coneAperture./max(coneAperture(:));





%% Red and green priamries

% Load in tables that give the shapes of the primary spectral radiance
% distributions, measured by the PR-650 
redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

% Change wavelength represnetation to be consistent with photoreceptors
Red_energy = SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl);
Green_energy = SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl);

% convert to quanta
Red_quantal = Red_energy.* (wvl.*1e-9)./(h*c);
Green_quantal = Green_energy.* (wvl.*1e-9)./(h*c);

% normalize to unit volume
Red_quantal = Red_quantal./sum(Red_quantal(:));
Green_quantal = Green_quantal./sum(Green_quantal(:));


LBkg = dot(L', normBkgRetIrrad);
MBkg = dot(M', normBkgRetIrrad);
SBkg = dot(S', normBkgRetIrrad);

LGreen = dot(L, Green_quantal);
LRed = dot(L, Red_quantal);
MGreen = dot(M, Green_quantal);
MRed = dot(M, Red_quantal);

% get thresholds

convertThresholdsToQuanta;

retIrradRedQuantaPerSecUm2 = retIrradRedQuantaPerSecM2 .* 1e-12;
retIrradGreenQuantaPerSecUm2 = retIrradGreenQuantaPerSecM2 .* 1e-12;

coneContrasts.red.dL_L = (LRed.*retIrradRedQuantaPerSecUm2)./(bkgScalar.*LBkg);
coneContrasts.red.dM_M = (MRed.*retIrradRedQuantaPerSecUm2)./(bkgScalar.*MBkg);

coneContrasts.green.dL_L = (LGreen.*retIrradGreenQuantaPerSecUm2)./(bkgScalar.*LBkg);
coneContrasts.green.dM_M = (MGreen.*retIrradGreenQuantaPerSecUm2)./(bkgScalar.*MBkg);

end

