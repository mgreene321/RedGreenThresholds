function [threshold543, threshold680, propL] = IdealObserverAnalysis2(coneData, stimLoc)

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
White = cal.rgb_spectra(whiteRGB(1),:,1) + cal.rgb_spectra(whiteRGB(2),:,2) + cal.rgb_spectra(whiteRGB(3),:,3);
[projectorRadianceWattsPerM2SrLambda, projectorRadianceS] = LumToRadiance(White', WlsToS(cal.wavelength_sampling'), projectorLum);
projectorRetIrradianceWattsPerM2Lambda = RadianceAndPupilAreaEyeLengthToRetIrradiance(projectorRadianceWattsPerM2SrLambda, projectorRadianceS, pupilAreaCm2,eyeLengthCm);
projectorRetIrradianceQuantaPerLambdaUm2Sec = projectorRetIrradianceWattsPerM2Lambda .* 1e-12 .* (cal.wavelength_sampling' .* 1e-9)/(h*c);
projectorRetIrradianceQuantaPerLambdaUm2Sec = SplineCmf(cal.wavelength_sampling', projectorRetIrradianceQuantaPerLambdaUm2Sec', SToWls(S));
for rr = 1:3
    T_quantalIsomerizations = photoreceptors.isomerizationAbsorptance(rr,:);
    ISareaUm2 = pi*(ISdiameter/2)^2;
    % isomerizations [1/sec/receptor] = inner segment area [um^2/receptor] * dot((quantal efficiency [1/quanta/lambda] * absorptance []), retinal irradiance [quanta/lambda/um^2/sec])
    projectorIsomerizationsSec(rr,:) = ISareaUm2*dot(T_quantalIsomerizations,projectorRetIrradianceQuantaPerLambdaUm2Sec);
    projectorFractionBleachedFromIsom(rr,:) = ComputePhotopigmentBleaching(projectorIsomerizationsSec(rr,:),'cones','isomerizations','Boynton');
end

photoreceptorsBleached = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptorsBleached.nomogram.lambdaMax = photoreceptors.nomogram.lambdaMax;

photoreceptorsBleached.fractionPigmentBleached.value = projectorFractionBleachedFromIsom;
photoreceptorsBleached = FillInPhotoreceptors(photoreceptorsBleached);

% At this point we have an image of the cone mosaic

imSize = [712 712];

%% Slice cone data

coneLocs = str2double(coneData(:, 1:2));
coneLabels = coneData(:,3);

% ignore missing cones

missingConeIdx = any(isnan(coneLocs),2) | any(coneLocs < 0,2) | ~(coneLabels == "L" | coneLabels == "M" | coneLabels == "S");

coneLocs = coneLocs(~missingConeIdx,:);
coneLabels = coneLabels(~missingConeIdx);

numCones = size(coneLocs,1);

%% Create cone aperture function

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
isoBkgPerLambda = nan(numCones,numWvls);

parfor lambda = 1:numWvls
    convBkg = projectorRetIrradianceQuantaPerLambdaUm2Sec(lambda).*pixelAreaUm2.*convNormBkgs(:,:,lambda); % each pixel gives quanta/s
    for cone = 1:numCones
        isoBkgPerLambda(cone,lambda) = photoreceptorsBleached.isomerizationAbsorptance(coneLabel2Idx(coneLabels(cone)),lambda) .* convBkg(coneLocs(cone,2), coneLocs(cone,1)) .* bkgDurSec;
        
%         switch coneLabels(cone)
%             case "L"
%                 lightFluxAtCone = convBkg(coneLocs(cone,2), coneLocs(cone,1));
%                 quantalEff = photoreceptorsBleached.isomerizationAbsorptance(1, lambda);
%                 %isoBkgPerLambda(cone,lambda) = quantalEff .* lightFluxAtCone .* bkgDurSec
%             case "M"
%                 quantalEff = photoreceptorsBleached.isomerizationAbsorptance(2, lambda);
%             case "S"
%                 quantalEff = photoreceptorsBleached.isomerizationAbsorptance(3,lambda);
%             otherwise
%                 quantalEff = coneProps(1).*photoreceptorsBleached.isomerizationAbsorptance(1,lambda)+...
%                     coneProps(2).*photoreceptorsBleached.isomerizationAbsorptance(2,lambda)+...
%                     coneProps(3).*photoreceptorsBleached.isomerizationAbsorptance(3,lambda)+...
%         end
    
    end
end

isoBkg = sum(isoBkgPerLambda,2);


%% Compute photon absorptions due to stimulus, assuming an intensity of 1;
greenFWHM = 26.2871; % nm
redFWHM = 30.6503; % nm
% 543
sigma = greenFWHM/(2*sqrt(2*log(2)));
normStimSpectrum = normpdf(wvl, 543, sigma);
normStimSpectrum = photoreceptorsBleached.preReceptoral.transmittance.*normStimSpectrum;
normStimSpectrum =normStimSpectrum./sum(normStimSpectrum(:));

% for each wavelengths
isoStimPerLambda543 = zeros(numCones,numWvls);

stimX = stimLoc(1);
stimY = stimLoc(2);

% compute stimulus isomerizations
parfor lambda = 1:numWvls
    
    shiftedStim = circshift(convNormStims(:,:,lambda), [stimY stimX] - [xc yc]);
    convStim = normStimSpectrum(lambda).*pixelAreaUm2.*shiftedStim.*stimDurSec;
    
    for cone = 1:numCones
        
        if shiftedStim(coneLocs(cone,2), coneLocs(cone,1)) >= 0.001
            
            
            isoStimPerLambda543(cone, lambda) = photoreceptorsBleached.quantalEfficiency.value(coneLabel2Idx(coneLabels(cone))).* photoreceptorsBleached.absorptance(coneLabel2Idx(coneLabels(cone)),lambda) .* convStim(coneLocs(cone,2), coneLocs(cone,1)) .* stimDurSec;
            %                 %isoBkg(cone,lambda) = isoBkgPerLambda;
            
            
            
            
        else
        end
        % end
    end
end



isoStim543 = sum(isoStimPerLambda543, 2);

%680
sigma = redFWHM/(2*sqrt(2*log(2)));
normStimSpectrum = normpdf(wvl, 680, sigma);
normStimSpectrum = photoreceptorsBleached.preReceptoral.transmittance.*normStimSpectrum;
normStimSpectrum =normStimSpectrum./sum(normStimSpectrum(:));

% for each wavelengths
isoStimPerLambda680 = zeros(numCones,numWvls);

stimX = stimLoc(1);
stimY = stimLoc(2);

parfor lambda = 1:numWvls
    
    shiftedStim = circshift(convNormStims(:,:,lambda), [stimY stimX] - [xc yc]);
    convStim = normStimSpectrum(lambda).*pixelAreaUm2.*shiftedStim.*stimDurSec;
    
    for cone = 1:numCones
        
        if shiftedStim(coneLocs(cone,2), coneLocs(cone,1)) >= 0.001
            
            
            isoStimPerLambda680(cone, lambda) = photoreceptorsBleached.quantalEfficiency.value(coneLabel2Idx(coneLabels(cone))).* photoreceptorsBleached.absorptance(coneLabel2Idx(coneLabels(cone)),lambda) .* convStim(coneLocs(cone,2), coneLocs(cone,1)) .* stimDurSec;
            %                 %isoBkg(cone,lambda) = isoBkgPerLambda;
        else
        end
        % end
    end
end

isoStim680 = sum(isoStimPerLambda680, 2);


%%
dPrimeThresh = 1.36;


%stimLoc = [304 244];
%stimLoc = [380 418];
f543 = @(stimIntensity) abs(dPrimeThresh - computeDprime(isoBkg,isoStim543, stimIntensity));
f680 = @(stimIntensity) abs(dPrimeThresh - computeDprime(isoBkg, isoStim680, stimIntensity));

options543 = optimset('TolX', 1e12, 'TolFun', 1e-3);
options680 = optimset('TolX', 1e14, 'TolFun', 1e-3);

threshold543 = fminsearch(f543, 1e17, options543);
threshold680 = fminsearch(f680, 1e19, options680);

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
