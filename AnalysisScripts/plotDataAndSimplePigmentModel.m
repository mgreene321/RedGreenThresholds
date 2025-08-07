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
ISdiameterUm = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)
stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
stimSideLengthPix = 21;
pixPerDeg = 560;
stimSideLengthDeg = stimSideLengthPix/pixPerDeg;
stimSideLengthUm = DegreesToRetinalMM(eyeLengthMm, stimSideLengthDeg) * 1e3;
stimAreaUm2 = stimSideLengthUm^2;
ISareaUm2 = pi.*(ISdiameterUm/2).^2;
%% Utility calculations
scanAngleDeg = 0.9;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
eyeLengthCm = eyeLengthMm*(10^-1);
scanAreaDegrees2 = scanAngleDeg^2;

%% Set up photoreceptors structure

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%photoreceptors.nomogram.lambdaMax = [555.5 530.3 420.7]';
photoreceptors.nomogram.lambdaMax = [563.4 530.3 420.7]';

photoreceptors = FillInPhotoreceptors(photoreceptors);

S = photoreceptors.nomogram.S;
wvl = SToWls(S);
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
%% Predicted red and green sensitivities

 greenFWHM = 26.2871; % nm
 redFWHM = 30.6503; % nm

 redSigma = redFWHM/(2*sqrt(2*log(2)));
 greenSigma = greenFWHM/(2*sqrt(2*log(2)));
% 
% Red = zeros(size(wvl));
% Red(wvl==680)=1;

% Green = zeros(size(wvl));
% Green(wvl==543)=1;

Red_energy = normpdf(wvl,680, redSigma);
Green_energy = normpdf(wvl, 543,greenSigma);

%   Red = Red./max(Red(:));
%   Green = Green./max(Green(:));

% convert to quanta
Red_quantal = Red_energy.* (wvl.*1e-9)./(h*c);
Green_quantal = Green_energy.* (wvl.*1e-9)./(h*c);

Red_quantal = Red_quantal./sum(Red_quantal(:));
Green_quantal = Green_quantal./sum(Green_quantal(:));

propLFine = linspace(0,1,1e3);

L = photoreceptorsBleached.isomerizationAbsorptance(1,:);
M = photoreceptorsBleached.isomerizationAbsorptance(2,:);

%L = photoreceptorsBleached.energyFundamentals(1,:);
%M = photoreceptorsBleached.energyFundamentals(2,:);
RedSensitivity = propLFine .* (dot(Red_quantal', L) - dot(Red_quantal', M)) + dot(Red_quantal', M);
GreenSensitivity = propLFine .* (dot(Green_quantal', L) - dot(Green_quantal', M)) + dot(Green_quantal', M);
%% Plot Data

propL = analyzedDataTable.L./(analyzedDataTable.L + analyzedDataTable.M);
redThresholdAUtimesRedPower = analyzedDataTable.RedThresholdPowerAU .* analyzedDataTable.RedPowerWattsAt680nm;
greenThresholdAUtimesGreenPower = analyzedDataTable.GreenThresholdPowerAU .* analyzedDataTable.GreenPowerWattsAt543nm;

for i = 1:size(analyzedDataTable,1)
    % For red stimuli
    redThresholdAU = analyzedDataTable.RedThresholdPowerAU(i);
    redPowerWattsAt680nm = analyzedDataTable.RedPowerWattsAt680nm(i);
    redThresholdRU = convertAOArbitraryToRealUnits(redThresholdAU,redPowerWattsAt680nm,Red_energy,S,scanAngleDeg, photoreceptorsBleached);
    
    incidentRedQuanta(i,:) = sum(redThresholdRU.rawCornIrradianceQuantaPerLambdaUm2SecIn) .* stimRenderTimeSec .* pupilAreaMm2 .* 1e6;
    transmittedRedQuanta(i,:) = sum(redThresholdRU.transmittedQuantaPerLambdaUm2Sec) .* stimRenderTimeSec;
 %   effectiveRedQuantaLMS(i,:) = sum(redThresholdRU.effectiveQuantaPerLambdaUm2Sec,2) * stimRenderTimeSec * ISareaUm2;
    
    % For green stimuli
    greenThresholdAU = analyzedDataTable.GreenThresholdPowerAU(i);
    greenPowerWattsAt543nm = analyzedDataTable.GreenPowerWattsAt543nm(i);
   greenThresholdRU = convertAOArbitraryToRealUnits(greenThresholdAU,greenPowerWattsAt543nm,Green_energy,S,scanAngleDeg, photoreceptorsBleached);
    
    incidentGreenQuanta(i,:) = sum(greenThresholdRU.rawCornIrradianceQuantaPerLambdaUm2SecIn) .* stimRenderTimeSec .* pupilAreaMm2 .* 1e6;
    transmittedGreenQuanta(i,:) = sum(greenThresholdRU.transmittedQuantaPerLambdaUm2Sec) .* stimRenderTimeSec;
%    effectiveGreenQuantaLMS(i,:) = sum(greenThresholdRU.effectiveQuantaPerLambdaUm2Sec,2) * stimRenderTimeSec * ISareaUm2;
    
end

figure; hold on
gscatter(propL, greenThresholdAUtimesGreenPower./redThresholdAUtimesRedPower, analyzedDataTable.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Green threshold (au) * Green Power / Red threshold (au) * RedPower')
ylim([0 0.05])

figure; hold on
gscatter(propL, incidentGreenQuanta./incidentRedQuanta, analyzedDataTable.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Incident green quanta / incident red quanta')
ylim([0 0.05])

figure; hold on
gscatter(propL, transmittedGreenQuanta./transmittedRedQuanta, analyzedDataTable.SubjectID, [], 'o');
plot(propLFine, RedSensitivity./GreenSensitivity, 'k--', 'LineWidth', 2, 'DisplayName', 'Predicted R/G Sensitivity')
xlabel('L/(L+M)')
ylabel('Red : Green Sensitivity');
title('Transmitted green quanta / Transmitted red quanta')
ylim([0 0.05])

tbl = table(analyzedDataTable.SubjectID, propL, redThresholdAUtimesRedPower, greenThresholdAUtimesGreenPower, incidentRedQuanta, incidentGreenQuanta, transmittedRedQuanta, transmittedGreenQuanta);
tbl.Properties.VariableNames = {'SubjectID', 'ProportionL', 'redThresholdAUtimesRedPower', 'greenThresholdAUtimesGreenPower', 'incidentRedQuantaPerReceptor', 'incidentGreenQuantaPerReceptor', 'transmittedRedQuantaPerReceptor', 'transmittedGreenQuantaPerReceptor'};



