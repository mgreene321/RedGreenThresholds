function [RedSensitivity, GreenSensitivity] = simplePigmentModel(params, propLdomain, wL)

if nargin < 3
    wL = 1;
end

%% Constants

loadConstants;
% h = 6.62607015e-34; % Planck's constant
% c = 2.99792458e8; % Speed of light
% pupilDiamMm = 6.5;
% pupilDiamCm = pupilDiamMm*(10^-1);
% eyeLengthMm = 16.7;
% pupilAreaCm2 = pi*((pupilDiamCm/2)^2);
% eyeLengthCm = eyeLengthMm*(10^-1);

%% Set up photoreceptors structure
photoreceptors = photoreceptorSensitivitiesFromScratch(params);

%photoreceptors = customConeSensitivities(params);

S = photoreceptors.nomogram.S;
wvl = SToWls(S);
%% Account for bleaching by projector background

projectorParams.cal = cal;
eyeParams.pupilAreaCm2 = pupilAreaCm2;
eyeParams.eyeLengthCm = eyeLengthCm;
%eyeParams.ISdiameterUm = coneDiamUm;
eyeParams.coneAperture = coneAperture;

pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;


% Background retinal power per pixel
projectorRetIrradianceQuantaPerLambdaPix2Sec = photoreceptors.projectorRetIrradianceQuantaPerSecUm2lambda * Um2PerPix;
normBkgRetIrrad = projectorRetIrradianceQuantaPerLambdaPix2Sec./sum(projectorRetIrradianceQuantaPerLambdaPix2Sec(:));

% L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
% M = photoreceptorsBleached.isomerizationAbsorptance(2,:);
% S = photoreceptorsBleached.isomerizationAbsorptance(3,:);

L = photoreceptors.isomerizationAbsorptanceBleached(1,:);
M = photoreceptors.isomerizationAbsorptanceBleached(2,:);
S = photoreceptors.isomerizationAbsorptanceBleached(3,:);

L = SplineCmf(SToWls(photoreceptors.nomogram.S), L, cal.wavelength_sampling',2);
M = SplineCmf(SToWls(photoreceptors.nomogram.S), M, cal.wavelength_sampling',2);
S = SplineCmf(SToWls(photoreceptors.nomogram.S), S, cal.wavelength_sampling',2);


%% Predicted red and green sensitivities

redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

% was splineCmf
%Red_energy = SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl);
%Green_energy = SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl);

Red_energy = redPrimaryTbl.power;
Green_energy = greenPrimaryTbl.power;
wvl = redPrimaryTbl.wavelength;

% convert to quanta
Red_quantal = Red_energy.* (wvl.*1e-9)./(h*c);
Green_quantal = Green_energy.* (wvl.*1e-9)./(h*c);

Red_quantal = Red_quantal./sum(Red_quantal(:));
Green_quantal = Green_quantal./sum(Green_quantal(:));



% L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
% M = photoreceptorsBleached.isomerizationAbsorptance(2,:);
% 
% L = photoreceptors.isomerizationAbsorptance(1,:);
% M = photoreceptors.isomerizationAbsorptance(2,:);

LRed = dot(Red_quantal', L);
LGreen = dot(Green_quantal', L);

MRed = dot(Red_quantal', M);
MGreen = dot(Green_quantal', M);

% 
% RedSensitivity = propLFine .* (dot(Red_quantal', L) - dot(Red_quantal', M)) + dot(Red_quantal', M);
% GreenSensitivity = propLFine .* (dot(Green_quantal', L) - dot(Green_quantal', M)) + dot(Green_quantal', M);

RedSensitivity = propLdomain.*wL.*LRed + (1-propLdomain).*MRed;
GreenSensitivity = propLdomain.*wL.*LGreen + (1-propLdomain).*MGreen;

end