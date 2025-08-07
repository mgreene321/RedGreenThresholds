root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
analyzedDataTable = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);

%% Constants
loadConstants;
 S = [380 1 401];
 wvl = SToWls(S);

redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

% was SplineCMF
Red_energy = transpose(SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl));
Green_energy = transpose(SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl));


% % convert to quanta
Red_quantal = Red_energy.* (wvl'.*1e-9)./(h*c);
Green_quantal = Green_energy.* (wvl'.*1e-9)./(h*c);


for i = 1:size(analyzedDataTable,1)
    % For red stimuli
    redThresholdAU = analyzedDataTable.RedThresholdPowerAU(i);
    redPowerWattsAt680nm = 0.4.*analyzedDataTable.RedPowerWattsAt680nm(i);
    redThresholdRU = convertAOArbitraryToRealUnits(redThresholdAU,redPowerWattsAt680nm,680, Red_energy,S,scanAngleDeg);
    retIrradRedQuantaPerSecM2(i,:) = sum(redThresholdRU.retIrradianceQuantaPerSecM2lambda);
    redQuantaPerSecLambdaIntoEye(i,:) = sum(redThresholdRU.quantaPerSecLambdaIntoEye);

    % For green stimuli
    greenThresholdAU = analyzedDataTable.GreenThresholdPowerAU(i);
    greenPowerWattsAt543nm = 0.4.*analyzedDataTable.GreenPowerWattsAt543nm(i);
    greenThresholdRU = convertAOArbitraryToRealUnits(greenThresholdAU,greenPowerWattsAt543nm, 543, Green_energy,S,scanAngleDeg);

    %incidentGreenQuantaMax(i,:) = max(g
    retIrradGreenQuantaPerSecM2(i,:) = sum(greenThresholdRU.retIrradianceQuantaPerSecM2lambda);
    greenQuantaPerSecLambdaIntoEye(i,:) = sum(greenThresholdRU.quantaPerSecLambdaIntoEye);
end