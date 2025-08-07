function [stim, bkg] = estimateIsomerizationsSingleLocation(coneApertures, stimLoc, thresholdAU, maxLaserPowerWatts, lambdaMax, relativeSpectrum, relativeSpectrumS)
% Isomerizations due to background

% Load in 'analyzedDataTable', which contains the variables, who's
% variables are described in AnalyzedDataTable.m
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
analyzedDataTable = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);

%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
pupilDiamMm = 6.5;
eyeLengthMm = 17;
%ISdiameterUm = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)
%stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
nominalStimDur = 2*(1/30);
stimSideLengthPix = 21;
pixPerDeg = 560;
stimSideLengthDeg = stimSideLengthPix/pixPerDeg;
stimSideLengthUm = DegreesToRetinalMM(eyeLengthMm, stimSideLengthDeg) * 1e3;
stimAreaUm2 = stimSideLengthUm^2;

%% Utility calculations
scanAngleDeg = 0.9;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);
scanAreaDegrees2 = scanAngleDeg^2;

% pixel size
ppd = 560;
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;


%% Stimulus parameters
fps = 30;
pixPerSec = 20e6;
stimDurFrames = 3;
stimDiamPix = 21;
stimDurSec = (1/pixPerSec) * stimDiamPix^2 * stimDurFrames;
bkgDurSec = (stimDurFrames-1)/fps;

%% Set up photoreceptors structure

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%photoreceptors.nomogram.lambdaMax = [555.5 530.3 420.7]';
%photoreceptors.nomogram.lambdaMax = [563.4 530.3 420.7]';

photoreceptors = FillInPhotoreceptors(photoreceptors);

Swvl = photoreceptors.nomogram.S;
wvl = SToWls(Swvl);
numWvls = numel(wvl);

%% Red and green priamries


% redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
% greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');
% 
% % was splineCmf
% Red_energy = SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl);
% Green_energy = SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl);
% 
% % convert to quanta
% Red_quantal = Red_energy.* (wvl.*1e-9)./(h*c);
% Green_quantal = Green_energy.* (wvl.*1e-9)./(h*c);
% 
% Red_quantal = Red_quantal./sum(Red_quantal(:));
% Green_quantal = Green_quantal./sum(Green_quantal(:));

thresholdRU = convertAOArbitraryToRealUnits(thresholdAU, maxLaserPowerWatts, lambdaMax, relativeSpectrum, relativeSpectrumS, scanAngleDeg);
retIrradQuantaPerSecM2 = sum(thresholdRU.retIrradianceQuantaPerSecM2lambda);
retIrradQuantaPerSecPix = retIrradQuantaPerSecM2 .* 1e-12 .* Um2PerPix;
stimScalar = retIrradQuantaPerSecPix;
%% PSF parameters

ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05;


%% Slice cone data

coneLocs = str2double(coneData(:, 1:2));
coneLabels = coneData(:,3);

% ignore missing cones

missingConeIdx = any(isnan(coneLocs),2) | any(coneLocs < 0,2);

coneLocs = coneLocs(~missingConeIdx,:);
coneLabels = coneLabels(~missingConeIdx);

numCones = size(coneLocs,1);

whiteRGB = [82 90 128];
%ND = -log10(20.39/170.2);
ND = 1;
projectorLum = 170.2 * 10^ND;

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
eyeParams.ISdiameterUm = coneDiamUm;

[photoreceptorsBleached, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);
projectorRetIrradianceQuantaPerLambdaPixSec = projectorRetIrradianceQuantaPerLambdaUm2Sec .* Um2PerPix;

L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
M = photoreceptorsBleached.isomerizationAbsorptance(2,:);
S = photoreceptorsBleached.isomerizationAbsorptance(3,:);

PSFfileName = fullfile(RedGreenThresholdsPath, ['PSFs_' num2str(S(1)) '_' num2str(S(2)) '_' num2str(S(3)) '.mat']);

if isfile(PSFfileName)
    PSF = importdata(PSFfileName);
else
    parfor lambda = 1:numWvls
        lambdaUm = wvl(lambda)/1000;
        tempPSF = generate_PSF(ap_field, psf_pupil, zernike_pupil, field_size, lambdaUm, diff_limited, defocus);
        
        tempPSF = tempPSF./sum(tempPSF(:)); % normalize volume to 1
        tempPSF = padarray(tempPSF, [100 100], 0, 'both');
        PSF(:,:,lambda) = tempPSF;
    end
    save(PSFfileName, 'PSF');
end

imSize = 712;
% normBkgs = ones(imSize-200, imSize-200, numel(wvl));
% normBkgs = padarray(normBkgs, [100 100], 0.5, 'both'); % 1 quanta/s/pixel
%normBkgs = normBkgs*Um2PerPix; % pixelAreaUm2 = area of pixel in um2 = 0.2807, so 1 quanta/s/pixel * pixelAreaUm2 [um^2/pixel] = 0.2807 quanta/s/um^2

normBkg = ones(imSize-200, imSize-200);
normBkg = padarray(normBkg, [100 100], 0.5, 'both'); % 1 quanta/s/pixel
%% Generate or load in bank of stimuli convolved with PSF and cone aperture function

% Stimulus intensity = 1 within boundary, 0 elsewhere 
stimTemplate = zeros(imSize, imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;

% Convolve 543 and 680 nm stimuli with PSF and aperture transfer function
convNormStim = fftshift(ifft2(fft2(stimTemplate).*fft2(PSF(:,:, wvl==lambdaMax))));


normBkgRetIrrad = projectorRetIrradianceQuantaPerLambdaPixSec./sum(projectorRetIrradianceQuantaPerLambdaPixSec(:));
bkgScalar = sum(projectorRetIrradianceQuantaPerLambdaPixSec(:));

% Cone sensitivities to the projector background
LBkg = dot(L', normBkgRetIrrad);
MBkg = dot(M', normBkgRetIrrad);
SBkg = dot(S', normBkgRetIrrad);


bkg_L = normBkg .* bkgScalar.*LBkg .* nominalStimDur;
bkg_M = normBkg .* bkgScalar.*MBkg .* nominalStimDur;
bkg_S = normBkg .* bkgScalar.*SBkg .* nominalStimDur;

% LGreen = dot(L, Green_quantal);
% LRed = dot(L, Red_quantal);
% MGreen = dot(M, Green_quantal);
% MRed = dot(M, Red_quantal);
% SGreen = dot(S,Green_quantal);
% SRed = dot(S,Red_quantal);

LStim = dot(L, relativeSpectrum);
MStim = dot(M, relativeSpectrum);
SStim = dot(S, relativespectrum);

% for each wavelengths

stimX = stimLoc(1);
stimY = stimLoc(2);

% compute stimulus isomerizations
xc = 357; yc = 357;
shiftedStim = circshift(convNormStim, [stimY stimX] - [xc yc]);

%%%%%

shiftedStim_L = shiftedStim .*stimScalar .* LStim.*nominalStimDur;
shiftedStim_M = shiftedStim .*stimScalar.* MStim.*nominalStimDur;
shiftedStim_S = shiftedStim .*stimScalar .* SStim.*nominalStimDur;

%% 
coneApertures_minPropL = coneApertures;
coneApertures_minPropL(:,:,1) = max(coneApertures(:,:,1), coneApertures(:,:,4));
coneApertures_minPropL(:,:,4) = [];

coneApertures_maxPropL = coneApertures;
coneApertures_maxPropL(:,:,2) = max(coneApertures(:,:,2), coneApertures(:,:,4));
coneApertures_maxPropL(:,:,4) = [];

% final background iso
Liso_bkg = bkg_L .* coneApertures(:,:,1);
Miso_bkg = bkg_M .* coneApertures(:,:,2);
Siso_bkg = bkg_S .* coneApertures(:,:,3);

Liso_bkg_minPropL = bkg_L .* coneApertures_minPropL(:,:,1);
Miso_bkg_minPropL = bkg_M .* coneApertures_minPropL(:,:,2);
Siso_bkg_minPropL = bkg_S .* coneApertures_minPropL(:,:,3);

Liso_bkg_maxPropL = bkg_L .* coneApertures_maxPropL(:,:,1);
Miso_bkg_maxPropL = bkg_M .* coneApertures_maxPropL(:,:,2);
Siso_bkg_maxPropL = bkg_S .* coneApertures_maxPropL(:,:,3);

% stim iso

Liso_Stim = shiftedStim_L .* coneApertures(:,:,1);
Miso_Stim = shiftedStim_M .* coneApertures(:,:,2);
Siso_Stim = shiftedStim_S .* coneApertures(:,:,3);

Liso_Stim_minPropL = shiftedStim_L .* coneApertures_minPropL(:,:,1);
Miso_Stim_minPropL = shiftedStim_M .* coneApertures_minPropL(:,:,2);
Siso_Stim_minPropL = shiftedStim_S .* coneApertures_minPropL(:,:,3);

Liso_Stim_maxPropL = shiftedStim_L .* coneApertures_maxPropL(:,:,1);
Miso_Stim_maxPropL = shiftedStim_M .* coneApertures_maxPropL(:,:,2);
Siso_Stim_maxPropL = shiftedStim_S .* coneApertures_maxPropL(:,:,3);

%%

iso_bkg = sum(Liso_bkg + Miso_bkg + Siso_bkg, 'all');
iso_bkg_minPropL = sum(Liso_bkg_minPropL + Miso_bkg_minPropL + Siso_bkg_minPropL, 'all');
iso_bkg_maxPropL = sum(Liso_bkg_maxPropL + Miso_bkg_maxPropL + Siso_bkg_maxPropL, 'all');

iso_Stim = sum(Liso_Stim + Miso_Stim + Siso_Stim, 'all');
iso_Stim_minPropL = sum(Liso_Stim_minPropL + Miso_Stim_minPropL + Siso_Stim_minPropL, 'all');
iso_Stim_maxPropL = sum(Liso_Stim_maxPropL + Miso_Stim_maxPropL + Siso_Stim_maxPropL, 'all');

bkg.iso_bkg = iso_bkg;
bkg.iso_bkg_minPropL = iso_bkg_minPropL;
bkg.iso_bkg_maxPropL = iso_bkg_maxPropL;

stim.iso_Stim = iso_Stim;
stim.iso_Stim_minPropL = iso_Stim_minPropL;
stim.iso_Stim_maxPropL = iso_Stim_maxPropL;

