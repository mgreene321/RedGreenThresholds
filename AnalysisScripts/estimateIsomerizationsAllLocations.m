root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);
Swvl = [380 1 401];
wvl = SToWls(Swvl);
redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
pupilDiamMm = 6.5;
eyeLengthMm = 17; %17;
ISdiameter = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)



% was splineCmf
Red_energy = SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl);
Green_energy = SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl);
% convert to quanta


Red_quantal = Red_energy.* (wvl.*1e-9)./(h*c);
Green_quantal = Green_energy.* (wvl.*1e-9)./(h*c);
Red_quantal = Red_quantal./sum(Red_quantal(:));
Green_quantal = Green_quantal./sum(Green_quantal(:));

for i = 1:size(ADT,1)
    subjectId = ADT.SubjectID(i);
    expFolder = ADT.Folder(i);
    stimLoc = [ADT.meanXc(i) ADT.meanYc(i)];
    stimLoc = round(stimLoc);
    redPower = ADT.RedPowerWattsAt680nm(i);
    greenPower = ADT.GreenPowerWattsAt543nm(i);
    
    redThresholdAU = 0.4*ADT.RedThresholdPowerAU(i);
    greenThresholdAU = 0.4*ADT.GreenThresholdPowerAU(i);
   coneAperturesFile = fullfile(RedGreenThresholdsPath,subjectId, expFolder, 'TransferFromOCT', 'coneApertures.mat');
   load(coneAperturesFile{1});
   
   [stimRed, bkg] = estimateIsomerizationsSingleLocation(coneApertures, stimLoc, redThresholdAU, redPower, 680, Red_quantal, Swvl);
   [stimGreen, ~] = estimateIsomerizationsSingleLocation(coneApertures, stimLoc, greenThresholdAU, greenPower, 543, Green_quantal, Swvl);
    
   bkg_L(i) = bkg.L;
   bkg_M(i) = bkg.M;
   bkg_S(i) = bkg.S;
   
   stimRed_L(i) = stimRed.L;
   stimRed_M(i) = stimRed.M;
   stimRed_S(i) = stimRed.S;
   
   stimGreen_L(i) = stimGreen.L;
   stimGreen_M(i) = stimGreen.M;
   stimGreen_S(i) = stimGreen.S;
   
   
end