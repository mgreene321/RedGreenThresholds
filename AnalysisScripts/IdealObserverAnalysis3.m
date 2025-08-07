function [threshold543, threshold680, propL] = IdealObserverAnalysis3(coneData, stimLoc, coneProps)

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
ISdiameter = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)
stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
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
pixelAreaUm2 = pixelSideLengthUm.^2;

%% stimulus parameters
fps = 30;
stimDurFrames = 3;
stimDiamPix = 21;
stimDurSec = 50e-9 * stimDiamPix^2 * stimDurFrames;
bkgDurSec = (stimDurFrames-1)/fps;

%% Set up photoreceptors structure

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%photoreceptors.nomogram.lambdaMax = [555.5 530.3 420.7]';
photoreceptors.nomogram.lambdaMax = [563.4 530.3 420.7]';

photoreceptors = FillInPhotoreceptors(photoreceptors);

S = photoreceptors.nomogram.S;
wvl = SToWls(S);
numWvls = numel(wvl);


%% PSF parameters

ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05;

%% Account for bleaching by projector background

whiteRGB = [82 90 128];
ND = -log10(20.39/170.2);
projectorLum = 170.2 * 10^ND;

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
eyeParams.ISdiameterUm = ISdiameterUm;

photoreceptorsBleached = projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, photoreceptors);

% At this point we have an image of the cone mosaic

imSize = [712 712];

%% Slice cone data

coneLocs = str2double(coneData(:, 1:2));
coneLabels = coneData(:,3);

% ignore missing cones

missingConeIdx = any(isnan(coneLocs),2) | any(coneLocs < 0,2);

coneLocs = coneLocs(~missingConeIdx,:);
coneLabels = coneLabels(~missingConeIdx);

numCones = size(coneLocs,1);

%% Create cone aperture function, whose volume = 1

coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = imSize;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma);


%% Generate or load in bank of PSFs

% generate psf for each wavelength, or load in file if it exists

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

%% Generate or load in bank of backgrounds convolved with PSF and cone aperture function
imSize = 712;
bkgTemplate = ones(imSize, imSize);

convNormBkgsFileName = fullfile(RedGreenThresholdsPath, 'convNormBkgs.mat');

if isfile(convNormBkgsFileName)
    convNormBkgs = importdata(convNormBkgsFileName);
    
else
    parfor lambda = 1:numWvls
        convNormBkgs(:,:,lambda) = ifft2(fft2(bkgTemplate).*fft2(PSF(:,:,lambda)).*fft2(coneAperture));
    end
    save(convNormBkgsFileName, 'convNormBkgs');
end


%% Generate or load in bank of stimuli convolved with PSF and cone aperture function
%stimSize = 21;
stimTemplate = zeros(imSize, imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;

convNormStimsFileName = fullfile(RedGreenThresholdsPath, 'convNormStims.mat');

if isfile(convNormStimsFileName)
    convNormStims = importdata(convNormStimsFileName);
else
    parfor lambda = 1:numWvls
        convNormStims(:,:,lambda) = ifft2(fft2(stimTemplate).*fft2(PSF(:,:,lambda)).*fft2(coneAperture));
    end
    save(convNormStimsFileName, 'convNormStims');
end


%% cone counting as sanity check

% as a sanity check count cones of each type
Lcounter = 0;
Mcounter = 0;
Scounter = 0;
xc = 357; yc = 357;
shiftedStim = circshift(convNormStims(:,:,find(wvl==543)), [stimLoc(2) stimLoc(1)] - [xc yc]);

for cone = 1:numCones
    
    %  if ~isnan(coneLocs(cone,1)) && ~isnan(coneLocs(cone,2)) && coneLocs(cone,1) > 0 && coneLocs(cone,2) > 0
    
    if shiftedStim(coneLocs(cone,2), coneLocs(cone,1)) >= 0.25
        
        if coneLabels(cone) == "L"
            Lcounter = Lcounter + 1;
            
        elseif coneLabels(cone) == "M"
            Mcounter = Mcounter + 1;
            
        elseif coneLabels(cone) == "S"
            Scounter = Scounter + 1;
            
        else % missing label
            % isoStimPerLambda(cone,lambda) = nan;
        end
        %isoBkg(cone,lambda) = isoBkgPerLambda;
        
    else
    end
    % end
    
end

propL = Lcounter/(Lcounter+Mcounter);
%% Compute photon absorptions due to background
%isoBkgPerLambda = nan(numCones,numWvls);

convBkgs = pixelAreaUm2.*convNormBkgs;
     % scale each page of convBkgs by the appropriate retinal irradiance
retIrrads = repelem(projectorRetIrradianceQuantaPerLambdaUm2Sec, prod(size(convNormBkgs, 1:2)))';
convBkgs = reshape(convBkgs(:).*retIrrads, size(convBkgs));

convBkgs_L = reshape(convBkgs(:) .* repelem(photoreceptorsBleached.isomerizationAbsorptance(1,:), prod(size(convNormBkgs, 1:2)))',  size(convNormBkgs));
convBkgs_M = reshape(convBkgs(:) .* repelem(photoreceptorsBleached.isomerizationAbsorptance(2,:), prod(size(convNormBkgs, 1:2)))',  size(convNormBkgs));
convBkgs_S = reshape(convBkgs(:) .* repelem(photoreceptorsBleached.isomerizationAbsorptance(3,:), prod(size(convNormBkgs, 1:2)))', size(convNormBkgs));

convBkgs_L = sum(convBkgs_L,3);
convBkgs_M = sum(convBkgs_M,3);
convBkgs_S = sum(convBkgs_S,3);


parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            isoBkg(cone,:) = convBkgs_L(coneLocs(cone,2), coneLocs(cone,1));
        case "M"
            isoBkg(cone,:) = convBkgs_M(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            isoBkg(cone,:) = convBkgs_S(coneLocs(cone,2), coneLocs(cone,1));
        otherwise
            isoBkg(cone,:) = coneProps(1).*convBkgs_L(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(2).*convBkgs_M(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(3).*convBkgs_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end

%% Compute photon absorptions due to stimulus, assuming an intensity of 1;

Lfactor = photoreceptorsBleached.quantalEfficiency.value(1).*photoreceptorsBleached.absorptance(1,:);
Mfactor = photoreceptorsBleached.quantalEfficiency.value(2).*photoreceptorsBleached.absorptance(2,:);
Sfactor = photoreceptorsBleached.quantalEfficiency.value(3).*photoreceptorsBleached.absorptance(3,:);



greenFWHM = 26.2871; % nm
redFWHM = 30.6503; % nm
% 543
sigma = greenFWHM/(2*sqrt(2*log(2)));
normStimSpectrum = normpdf(wvl, 543, sigma);
normStimSpectrum = photoreceptorsBleached.preReceptoral.transmittance.*normStimSpectrum';
normStimSpectrum =normStimSpectrum./sum(normStimSpectrum(:));

% for each wavelengths

stimX = stimLoc(1);
stimY = stimLoc(2);

% compute stimulus isomerizations

shiftedStim = circshift(convNormStims, [stimY stimX] - [xc yc]);

%%%%%
shiftedStim = pixelAreaUm2.*shiftedStim;
     % scale each page of convBkgs by the appropriate retinal irradiance
retIrrads = repelem(normStimSpectrum, prod(size(shiftedStim, 1:2)))';
shiftedStim = reshape(shiftedStim(:).*retIrrads, size(shiftedStim));

shiftedStim_L = reshape(shiftedStim(:) .* repelem(Lfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_M = reshape(shiftedStim(:) .* repelem(Mfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_S = reshape(shiftedStim(:) .* repelem(Sfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_L = sum(shiftedStim_L,3);
shiftedStim_M = sum(shiftedStim_M,3);
shiftedStim_S = sum(shiftedStim_S,3);


parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            isoStim543(cone,:) = shiftedStim_L(coneLocs(cone,2), coneLocs(cone,1));
        case "M"
            isoStim543(cone,:) = shiftedStim_M(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            isoStim543(cone,:) = shiftedStim_S(coneLocs(cone,2), coneLocs(cone,1));
        otherwise
            isoStim543(cone,:) = coneProps(1).*shiftedStim_L(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(2).*shiftedStim_M(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(3).*shiftedStim_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end



%680
sigma = redFWHM/(2*sqrt(2*log(2)));
normStimSpectrum = normpdf(wvl, 680, sigma);
normStimSpectrum = photoreceptorsBleached.preReceptoral.transmittance.*normStimSpectrum';
normStimSpectrum =normStimSpectrum./sum(normStimSpectrum(:));


shiftedStim = circshift(convNormStims, [stimY stimX] - [xc yc]);

%%%%%
shiftedStim = pixelAreaUm2.*shiftedStim;
     % scale each page of convBkgs by the appropriate retinal irradiance
retIrrads = repelem(normStimSpectrum, prod(size(shiftedStim, 1:2)))';
shiftedStim = reshape(shiftedStim(:).*retIrrads, size(shiftedStim));

shiftedStim_L = reshape(shiftedStim(:) .* repelem(Lfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_M = reshape(shiftedStim(:) .* repelem(Mfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_S = reshape(shiftedStim(:) .* repelem(Sfactor, prod(size(shiftedStim, 1:2)))',  size(shiftedStim));
shiftedStim_L = sum(shiftedStim_L,3);
shiftedStim_M = sum(shiftedStim_M,3);
shiftedStim_S = sum(shiftedStim_S,3);


parfor cone = 1:numCones
    switch coneLabels(cone)
        case "L"
            isoStim680(cone,:) = shiftedStim_L(coneLocs(cone,2), coneLocs(cone,1));
        case "M"
            isoStim680(cone,:) = shiftedStim_M(coneLocs(cone,2), coneLocs(cone,1));
        case "S"
            isoStim680(cone,:) = shiftedStim_S(coneLocs(cone,2), coneLocs(cone,1));
        otherwise
            isoStim680(cone,:) = coneProps(1).*shiftedStim_L(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(2).*shiftedStim_M(coneLocs(cone,2), coneLocs(cone,1)) + coneProps(3).*shiftedStim_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end

%%%%%



%%
dPrimeThresh = 1.36;


%stimLoc = [304 244];
%stimLoc = [380 418];
f543 = @(stimIntensity) abs(dPrimeThresh - computeDprime(isoBkg,isoStim543, stimIntensity));
f680 = @(stimIntensity) abs(dPrimeThresh - computeDprime(isoBkg, isoStim680, stimIntensity));

options543 = optimset('TolX', 1e3, 'TolFun', 1e-3);
options680 = optimset('TolX', 1e3, 'TolFun', 1e-3); %'Display','iter','PlotFcns',@optimplotfval);

threshold543 = fminsearch(f543, 1e4, options543);
threshold680 = fminsearch(f680, 1e5, options680);

function dprime = computeDprime(backgroundIsomerizations, stimulusIsomerizations, stimIntensity)

stimulusIsomerizations = stimIntensity.*stimulusIsomerizations + backgroundIsomerizations;

alpha = backgroundIsomerizations;
beta = stimulusIsomerizations;

dprime = sum((beta-alpha).*log(beta./alpha))./ sqrt(0.5.*(sum((beta+alpha).*(log(beta./alpha).^2))));

function coneIdx = coneLabel2Idx(coneLabel)

if coneLabel == "L"
    coneIdx = 1;
elseif coneLabel == "M"
    coneIdx = 2;
elseif coneLabel == "S"
    coneIdx = 3;
else
end
