%% King-Smith Carden
% Eye length

h = 6.626070150000000e-34;
c =  2.99792458e8;

eyeLengthM = 17e-3;

%Pupil size

pupilDiamMm_KSC = 3;
pupilDiamM_KSC = pupilDiamMm_KSC*(10^-3);
pupilAreaMm2_KSC = pi*(pupilDiamMm_KSC/2)^2;
pupilAreaM2_KSC = pi*(pupilDiamM_KSC/2)^2;

% Degrees to meter conversion

MperDeg = DegreesToRetinalMM(1, 17)*1e-3;
degPerM = 1/MperDeg;
deg2PerM2 = degPerM^2;

load T_xyz1931.mat;
S = [380 1 401];
wvl = SToWls(S);
T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:), S);


% Threshold at 543 nm

 t543_quantaPerSecDeg2_KSC = 10^9.25;
  
 
 % convert to trolands
 
conversionFactor = 683*h*c*((180/pi)^2) * 1e6;

t543_trolands_KSC = conversionFactor*t543_quantaPerSecDeg2_KSC*T_vLambda(wvl==543)/(543e-9);

% Convert to luminance 

t543_CdPerM2_KSC = 683 * T_vLambda(wvl==543)* t543_WattsPerM2Sr_KSC;
t543_trolands_KSC = t543_CdPerM2_KSC*pupilAreaMm2_KSC;

%% Our data

%Pupil size

pupilDiamMm = 6.5;
pupilDiamM = pupilDiamMm*(10^-3);
pupilAreaMm2 = pi*(pupilDiamMm/2)^2;
pupilAreaM2 = pi*(pupilDiamM/2)^2;

% Scan
scanAngleDeg = 0.9;

% Green SPD

greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');

Green_energy = transpose(SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl));
Green_quantal = Green_energy.* (wvl'.*1e-9)./(h*c);

Green_energy = Green_energy/sum(Green_energy(:));
Green_quantal = Green_quantal/sum(Green_quantal(:));

% Mean 543 nm power 
meanGreenPower = 8.785e-8; % will
meanGreenPower = 7.677e-8; % max

% Mean 543 nm threshold in au
meanGreenThresholdAU = 0.113; %will
meanGreenThresholdAU = 0.312; %max

t543_realUnits = convertAOArbitraryToRealUnits(meanGreenThresholdAU, meanGreenPower, 543,Green_quantal, S, scanAngleDeg);
t543_CdPerM2 = 683 * dot(T_vLambda, t543_realUnits.radianceWattsPerM2Srlambda);
t543_trolands= t543_CdPerM2*pupilAreaMm2;



