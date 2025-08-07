function figureIdealObserver543

fontSize = 16;
lineWidth = 2;
clear parula;
% pinkColorMap = parula;
% 
% pinkColorMapHSV = rgb2hsv(pinkColorMap);
% 
% pinkRGB = @(x) interp1(linspace(0,1,256), pinkColorMap, x);

numColorMapRows = 256;
pixVals = 1:256;

hsvMatrix_Green1 = zeros(numColorMapRows,3);
hsvMatrix_Green1(:,1) = 120/360; % green hue
hsvMatrix_Green1(:,2) =1;%pixVals.^-0.07; % vary saturation
hsvMatrix_Green1(:,3) = 0.7.*(pixVals.^(1/8) - 1) + 0.3.*linspace(0,1, numColorMapRows);%pinkColorMapHSV(:,3);%linspace(0.75,1, 256); %linspace(0.25,1,256); % keep "value" constant
rgbMatrix_Green1 = hsv2rgb(hsvMatrix_Green1);

greenRGB1 = @(x) interp1(linspace(0,1,256), rgbMatrix_Green1, x);

hsvMatrix_Green2 = zeros(numColorMapRows,3);
hsvMatrix_Green2(:,1) = 120/360; % red hue
hsvMatrix_Green2(:,2) = 0.9.*(pixVals/256)+0.1;%pixVals.^-0.07; % vary saturation
hsvMatrix_Green2(:,3) = 0.7.*(pixVals.^(1/8) - 1) + 0.3.*linspace(0,1, numColorMapRows);%pinkColorMapHSV(:,3);%linspace(0.75,1, 256); %linspace(0.25,1,256); % keep "value" constant
rgbMatrix_Green2 = hsv2rgb(hsvMatrix_Green2);

greenRGB2 = @(x) interp1(linspace(0,1,256), rgbMatrix_Green2, x);


%% cropping
% crop everyothing to be a 9 x 9 arcmin rectangle around the stimulus
ppd = 560;
rectSize = 85;

cropRect = [357-fix(rectSize/2),  357-fix(rectSize/2), rectSize, rectSize];

%%
% Ideal obsever figure
addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\AnalysisScripts'));

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
imSize = 712;
S = [380 1 401];
wvl = SToWls(S);
PSFfileName = fullfile(RedGreenThresholdsPath, ['PSFs_' num2str(S(1)) '_' num2str(S(2)) '_' num2str(S(3)) '.mat']);
PSF = importdata(PSFfileName);

% load in demo aperture
coneAperturesFile = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\5_21_2024\TransferFromOCT\coneApertures.mat';
load(coneAperturesFile);
coneAperturesNorm = coneApertures./max(coneApertures(:));
%load in mean frame
meanFrameFile = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\5_21_2024\TransferFromOCT\meanFrame.mat';
load(meanFrameFile);

% load in cone data
dataFile = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\10001R\5_21_2024\TransferFromOCT\dailyData.mat';
load(dataFile);
coneData = dailyData;
coneLocs = str2double(coneData(:, 1:2));
coneLabels = coneData(:,3);

% ignore missing cones
missingLabel = ~(coneLabels == "L" | coneLabels == "M" | coneLabels == "S");
missingConeIdx = any(isnan(coneLocs),2) | any(coneLocs < 0,2) | missingLabel;
coneLocs = coneLocs(~missingConeIdx,:);
coneLabels = coneLabels(~missingConeIdx);

numCones = size(coneLocs,1);


%% PSF
% 
% figure;
% [X,Y] = meshgrid(1:712);
%imagesc(imcrop(PSF(:,:,wvl==680), cropRect)), colormap gray
% surf(imcrop(X, cropRect), imcrop(Y, cropRect), imcrop(PSF(:,:, wvl==680), cropRect), 'EdgeColor', [0.5 0.5 0.5]); colormap gray
% xlim([cropRect(1), cropRect(1) + cropRect(3)]);
% ylim([cropRect(2), cropRect(2) + cropRect(4)]);
% set(gca, 'XColor', 'k', 'YColor', 'k', 'LineWidth', lineWidth, 'FontSize', fontSize);
% set(gca, 'color', 'none'), grid off;

figure;
imagesc(imcrop(PSF(:,:,wvl==543)./(max(max(PSF(:,:,(wvl==543))))), cropRect)), colormap(rgbMatrix_Green1), axis equal off



%% Light delivery contour
stimDiamPix = 21;

stimTemplate = zeros(imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;

stim543_blurred = fftshift(ifft2(fft2(stimTemplate).*fft2(PSF(:,:,wvl==543))));
lightWithinContours = [0.1 0.5 0.9];
levels = integratedLightContour(stim543_blurred, lightWithinContours);

stim543_blurred_crop = imcrop(stim543_blurred, cropRect);
C = contourc(stim543_blurred_crop./max(stim543_blurred_crop(:)), levels);
[xx, yy, ~, ~] = getContourCoordinates(C);

%figure, imagesc(imcrop(meanFrame, cropRect)), colormap gray, hold on;
redHue = 0;
% plot(xx{1}, yy{1}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 0.33]));
% plot(xx{2}, yy{2}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 0.67]));
% plot(xx{3}, yy{3}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 1]));

%% Light delivery on mosaic with cones labeled

figure, imagesc(imcrop(coneAperturesNorm(:,:,1:3), cropRect)); hold on; alpha(0.6); axis equal;

% plot(xx{1}, yy{1}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 0.33]));
% plot(xx{2}, yy{2}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 0.67]));
% plot(xx{3}, yy{3}, 'LineWidth', 3, 'Color', hsv2rgb([redHue 1 1]));
% % 
plot(xx{1}, yy{1}, 'LineWidth', 4, 'Color', greenRGB1(0.1));
plot(xx{2}, yy{2}, 'LineWidth', 4, 'Color', greenRGB1(0.6));
plot(xx{3}, yy{3}, 'LineWidth', 4, 'Color', greenRGB1(1));



%% Cone isomerizations

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
photoreceptors = FillInPhotoreceptors(photoreceptors);

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

%% isomerizations due to background
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



%% iso due to stimulus

LGreen = dot(L, Green_quantal);
LRed = dot(L, Red_quantal);
MGreen = dot(M, Green_quantal);
MRed = dot(M, Red_quantal);
SGreen = dot(S,Green_quantal);
SRed = dot(S,Red_quantal);

stimTemplate = zeros(imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;
stim543_centered = ifft2(fft2(stimTemplate).*fft2(PSF(:,:,wvl==543)).*fft2(coneAperture));

stimX = 357; stimY = 357;
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
%             isoStim680(cone,:) = coneProps(1).*shiftedStim680_L(coneLocs(cone,2), coneLocs(cone,1)) + ...
%                 coneProps(2).*shiftedStim680_M(coneLocs(cone,2), coneLocs(cone,1)) +...
%                 coneProps(3).*shiftedStim680_S(coneLocs(cone,2), coneLocs(cone,1));
    end
end

total543Isom = (10^9.7).*stim543Isom + bkgIsom;
isoms = zeros(712);

%cone = Circle(3.5);
for i = 1:numCones
    x = coneLocs(i,1); y = coneLocs(i,2);
    %isoms((y-3):(y+3),(x-3):(x+3)) = total680Isom(i).*cone;
    isoms(y,x) = total543Isom(i);
end

isoms = imcrop(isoms, cropRect);
isoms = isoms./max(isoms(:));

[b,a] = find(isoms);

canvas = ones(size(isoms));
figure; imagesc(canvas), colormap gray, hold on; axis equal;

for i =1:numel(a)
    %plot(a(i), b(i), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', isoms(b(i), a(i)).*ones(1,3), 'MarkerSize', 20)
    plot(a(i), b(i), 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', greenRGB2(isoms(b(i), a(i))), 'MarkerSize', 20)
end
 %colormap(rgbMatrix_Green2(6,:))
colormap([0 0 0]);
function pixVal = integratedLightContour(H, A)

 %normalize heatmap:
 
 H = H./max(H(:));
 
 h = max(max(abs(diff(H))));
 
 %pixelValues = 0:h:1;
 %pixelValues = unique(H);
 pixelValues = 0:0.01:1;
 
 if pixelValues(end)~=1
     pixelValues(end+1) = 1;
 else
 end
 
 
for i = 1:length(pixelValues)
    integratedLight(i) = sum(H(H>=pixelValues(i)))./sum(H(:));
    
end

% make sure sample points unique
[~, ia, ~] = intersect(integratedLight, unique(integratedLight));

I = @(x) interp1(integratedLight(ia), pixelValues(ia), x, 'linear', 'extrap');

if I(A) < 0
    pixVal = 0;
elseif I(A) > 1
    pixVal = 1;
else
    pixVal = I(A);
end
