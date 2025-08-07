
%% Physical constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light

%% Biological constants

% eye dimensions
pupilDiamMm = 6.5;
pupilDiamCm = pupilDiamMm*(10^-1);
eyeLengthMm = 16.7;

% cone aperture
coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = 11;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma);
%% Utility calculations
ppd = 560;
scanAngleDeg = 512/ppd;
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pi*((pupilDiamCm/2)^2);
eyeLengthCm = eyeLengthMm*(10^-1);

% projector parameters
whiteRGB = [82 90 128];
ND = 1;
projectorLum = (170.2 * 10^ND);
cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');

