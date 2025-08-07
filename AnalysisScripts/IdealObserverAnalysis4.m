function [threshold543, threshold680] = IdealObserverAnalysis4(coneData, stimLoc, coneProps)

% Get path to red green thresholds data
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);

%% Physical constants

h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
%% Eye constants

pupilDiamMm = 6.5;
eyeLengthMm = 16.7;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);

%% Resolution of imaging system
ppd = 560; % pixels per degree
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;
imSize = 712; % all the images generated below will be 712 x 712 pixels
%ISdiameterUm = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)

%% Stimulus parameters

stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
nominalStimDur = 2*(1/30); % time the stimulus is on based on the number of frames it is presented
stimDiamPix = 21;


%% Photoreceptors structure: L-cone lambda max is hard-coded here

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%photoreceptors.nomogram.lambdaMax = [555.5 530.3 420.7]';
%photoreceptors.nomogram.lambdaMax = [563.4 530.3 420.7]';
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Wavelength domain from photoreceptors structure
Swvl = photoreceptors.nomogram.S;
wvl = SToWls(Swvl);
numWvls = numel(wvl);

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

%% PSF parameters

ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05; 

%% Cone data

coneLocs = str2double(coneData(:, 1:2));
coneLabels = coneData(:,3);

% ignore missing cones
missingLabel = ~(coneLabels == "L" | coneLabels == "M" | coneLabels == "S");
missingConeIdx = any(isnan(coneLocs),2) | any(coneLocs < 0,2) | missingLabel;
coneLocs = coneLocs(~missingConeIdx,:);
coneLabels = coneLabels(~missingConeIdx);

numCones = size(coneLocs,1);

%% Create cone aperture function, whose volume = 1

coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 712;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma); % Gaussian cone aperture
coneAperture = coneAperture./max(coneAperture(:));

%% Account for bleaching by projector background

whiteRGB = [82 90 128];
%ND = -log10(20.39/170.2); % ND based on PR-650 measurement
ND = 1; % nominal ND
projectorLum = 170.2 * 10^ND;

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
eyeParams.coneAperture = fspecial('gaussian', 11, gaussianSigma); 

% Get cone quantum efficiencies at steady state bleached level, and
% projector retinal irradiance
[photoreceptorsBleached, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = ...
    projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);

% Background retinal power per pixel
projectorRetIrradianceQuantaPerLambdaPix2Sec = projectorRetIrradianceQuantaPerLambdaUm2Sec * Um2PerPix;
normBkgRetIrrad = projectorRetIrradianceQuantaPerLambdaPix2Sec./sum(projectorRetIrradianceQuantaPerLambdaPix2Sec(:));
bkgScalar = sum(projectorRetIrradianceQuantaPerLambdaPix2Sec(:));

L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
M = photoreceptorsBleached.isomerizationAbsorptance(2,:);
S = photoreceptorsBleached.isomerizationAbsorptance(3,:);

%% Generate or load in bank of PSFs

% generate psf for each wavelength, or load in file if it exists
PSFfileName = fullfile(RedGreenThresholdsPath, ['PSFs_' num2str(Swvl(1)) '_' num2str(Swvl(2)) '_' num2str(Swvl(3)) '.mat']);

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

%% ISOMERIZATIONS DUE TO BACKGROUND

% Generate background template
projectorPPD = 35;
patchSizeDeg = 40/projectorPPD;
patchSizeAOPixels = patchSizeDeg*ppd;

imSizeMinusPatchSize = imSize - patchSizeAOPixels;
normBkg = ones(patchSizeAOPixels);
normBkg = padarray(normBkg, [imSizeMinusPatchSize/2 imSizeMinusPatchSize/2], 0.5, 'both'); % 1 quanta/s/pixel

% Scale background by appropriate intensity
bkg = bkgScalar.*normBkg;
 % now each pixel gives the actual number of quanta/s

% Convolve background with cone aperture 
bkg = fftshift(ifft2(fft2(bkg).*fft2(coneAperture)));
 % now each pixel gives the quanta/s collected by a cone centered at that
 % pixel
 
% Compute cone sensitivities to background light
% these are the sensitivities of each cone to a hypothetical quanta of the background light 
LBkg = dot(L', normBkgRetIrrad);
MBkg = dot(M', normBkgRetIrrad);
SBkg = dot(S', normBkgRetIrrad);


% isomerizations/stimulus = quanta/s x isomerizations/quantum x s/stimulus
bkg_Lisom = bkg .*LBkg .* nominalStimDur;
bkg_Misom = bkg .*MBkg .* nominalStimDur;
bkg_Sisom = bkg .*SBkg .* nominalStimDur;

% Account for projector duration , and multiply by receptor area

parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            bkgIsom(cone,:) = bkg_Lisom(coneLocs(cone,2), coneLocs(cone,1)) ;
        case "M"
            bkgIsom(cone,:) = bkg_Misom(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            bkgIsom(cone,:) = bkg_Sisom(coneLocs(cone,2), coneLocs(cone,1));
%         otherwise
%             isoBkg(cone,:) = coneProps(1).*bkg_L(coneLocs(cone,2), coneLocs(cone,1))+...
%                 coneProps(2).*bkg_M(coneLocs(cone,2), coneLocs(cone,1)) + ...
%                 coneProps(3).*bkg_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end


%% ISOMERIZATIONS DUE TO STIMULUS
% Generate stimulus template
% Stimulus intensity = 1 within boundary, 0 elsewhere 
stimTemplate = zeros(imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;


stim543_centered = ifft2(fft2(stimTemplate).*fft2(PSF(:,:, wvl==543)).*fft2(coneAperture));
stim680_centered = ifft2(fft2(stimTemplate).*fft2(PSF(:,:,wvl==680)).*fft2(coneAperture));

%% Compute photon absorptions due to background
%isoBkgPerLambda = nan(numCones,numWvls);

% each pixel of convBkgs gives how much light would be abosrbed by a cone
% centered at that pixel, assuming a stimulus 

% Cone sensitivities to the projector background


%% Compute photon absorptions due to stimulus, assuming an intensity of 1;

LGreen = dot(L, Green_quantal);
LRed = dot(L, Red_quantal);
MGreen = dot(M, Green_quantal);
MRed = dot(M, Red_quantal);
SGreen = dot(S,Green_quantal);
SRed = dot(S,Red_quantal);

% for each wavelengths

stimX = stimLoc(1);
stimY = stimLoc(2);

% compute stimulus isomerizations
xc = 357; yc = 357;
stim543 = circshift(stim543_centered, [stimY stimX] - [xc yc]);

%%%%%

stim543_Lisom = stim543 .* LGreen.*stimRenderTimeSec;%.* coneAreaUm2;
stim543_Misom = stim543 .* MGreen.*stimRenderTimeSec;%.* coneAreaUm2;
stim543_Sisom = stim543 .* SGreen.*stimRenderTimeSec;%.* coneAreaUm2;

parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            stim543Isom(cone,:) = stim543_Lisom(coneLocs(cone,2), coneLocs(cone,1));
        case "M"
            stim543Isom(cone,:) = stim543_Misom(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            stim543Isom(cone,:) = stim543_Sisom(coneLocs(cone,2), coneLocs(cone,1));
%         otherwise
%             isoStim543(cone,:) = coneProps(1).*shiftedStim543_L(coneLocs(cone,2), coneLocs(cone,1)) +...
%                 coneProps(2).*shiftedStim543_M(coneLocs(cone,2), coneLocs(cone,1)) +...
%                 coneProps(3).*shiftedStim543_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end

%680

stim680 = circshift(stim680_centered, [stimY stimX] - [xc yc]);

%%%%%

stim680_Lisom = stim680 .* LRed.*stimRenderTimeSec;%.* coneAreaUm2;
stim680_Misom = stim680 .* MRed.*stimRenderTimeSec;%.* coneAreaUm2;
stim680_Sisom = stim680 .* SRed.*stimRenderTimeSec;%.* coneAreaUm2;


parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            stim680Isom(cone,:) = stim680_Lisom(coneLocs(cone,2), coneLocs(cone,1));
        case "M"
            stim680Isom(cone,:) = stim680_Misom(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            stim680Isom(cone,:) = stim680_Sisom(coneLocs(cone,2), coneLocs(cone,1));
%         otherwise
%             isoStim680(cone,:) = coneProps(1).*shiftedStim680_L(coneLocs(cone,2), coneLocs(cone,1)) + ...
%                 coneProps(2).*shiftedStim680_M(coneLocs(cone,2), coneLocs(cone,1)) +...
%                 coneProps(3).*shiftedStim680_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end

%%%%%
%%
dPrimeThresh = 1.36;

f543 = @(logStimIntensity) abs(dPrimeThresh - computeDprime(bkgIsom,stim543Isom, logStimIntensity));
f680 = @(logStimIntensity) abs(dPrimeThresh - computeDprime(bkgIsom, stim680Isom, logStimIntensity));

options = optimoptions('fmincon', 'StepTolerance', 1e-12);
threshold543 =fmincon(f543, 5, [], [], [], [], 0, Inf, [], options); %3.5
threshold680 =fmincon(f680,7, [], [], [], [], 0, Inf, [], options); %5.5

function dprime = computeDprime(backgroundIsomerizations, stimulusIsomerizations, logStimIntensity)

stimulusIsomerizations = (10^logStimIntensity).*stimulusIsomerizations + backgroundIsomerizations;

alpha = backgroundIsomerizations;
beta = stimulusIsomerizations;

dprime = sum((beta-alpha).*log(beta./alpha))./ sqrt(0.5.*(sum((beta+alpha).*(log(beta./alpha).^2))));

