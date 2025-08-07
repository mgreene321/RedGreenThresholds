function b = figureRawSensitivityAndVLambda

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);

%% Constants
loadConstants;
% h = 6.62607015e-34; % Planck's constant
% c = 2.99792458e8; % Speed of light
% pupilDiamMm = 6.5;
% pupilDiamCm = pupilDiamMm*(10^-1);
% eyeLengthMm = 16.7;
% ISdiameterUm = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)
% stimRenderTimeSec = 50e-9 * 21^2 * 3; % sec/pixel * number of pixels * number of frames
% stimSideLengthPix = 21;
% pixPerDeg = 560;
% stimSideLengthDeg = stimSideLengthPix/pixPerDeg;
% stimSideLengthUm = DegreesToRetinalMM(eyeLengthMm, stimSideLengthDeg) * 1e3;
% stimAreaUm2 = stimSideLengthUm^2;
% ISareaUm2 = pi.*(ISdiameterUm/2).^2;
% 
% ppd =560;
% coneSpacingPixels = 11;
% coneDiamPixels = coneSpacingPixels;
% coneDiamDeg = coneDiamPixels./ppd;
% coneDiamUm = DegreesToRetinalMM(coneDiamDeg, eyeLengthMm)*1e3;
% coneAreaUm2 = pi*(coneDiamUm/2)^2;
% %% Utility calculations
% scanAngleDeg = 0.9;
% pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
% pupilAreaCm2 = pi*((pupilDiamCm/2)^2);
% eyeLengthCm = eyeLengthMm*(10^-1);

axFontSize = 12;
labelFontSize = 14;
defaultLineWidth = 2;
defaultColor = [0 0 0]; %ie black
tickDir = 'out';
defaultMarker = 'o';
markerLineWidth = 1.5;
markerSize = 8;
tickLength = [0.025 0.025];
subjColors = [0 0 0.35; 1 0.4 0];
seed =1;
rng(seed);
lineStyles = {':', '-', '--'};


scanAngleDeg = 0.91;
scanAngleDegrees2 = scanAngleDeg^2;

% % to account for bleaching
% % Make cone aperture
% coneSpacingPixels = 11;
% coneDiamPixels = coneSpacingPixels;
% coneAperturePixels = coneDiamPixels.*.5;
% gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
% windowSize = 11;
% coneAperture = fspecial('gaussian', windowSize, gaussianSigma);
% 
% whiteRGB = [82 90 128];
% %ND = -log10(20.39/170.2);
% ND = 1;
% projectorLum = (170.2 * 10^ND);
% 
% cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
% projectorParams.cal = cal;
% eyeParams.pupilAreaCm2 = pupilAreaCm2;
% eyeParams.eyeLengthCm = eyeLengthCm;
% %eyeParams.ISdiameterUm = coneDiamUm;
% eyeParams.coneAperture = coneAperture;
% 
% ppd = 560; % pixels per degree
% pixelSideLengthDeg = 1/ppd;
% pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
% pixelSideLengthUm = pixelSideLengthMm*1e3;
% Um2PerPix= pixelSideLengthUm.^2;



figure; hold on
set(gcf, 'PaperPositionMode', 'auto', 'Renderer', 'painters', 'Color', 'w')
set(gca, 'Position', [0.225 0.2 0.8 0.8], 'FontSize', axFontSize, 'XColor', defaultColor, 'YColor', defaultColor, 'LineWidth', defaultLineWidth, 'TickDir', tickDir, 'box', 'off');
set(gca, 'TickLength', tickLength)

load T_CIE_Y2;
wls = SToWls(S_CIE_Y2);

%Convert V_lambda to quantal units
V_Lambda_quantal = EnergyToQuanta(wls, T_CIE_Y2');
V_Lambda_quantal = QuantaToEnergy(wls, T_CIE_Y2');
%get thresholds in quanta/s/deg2

convertThresholdsToQuanta;
tbl.SubjectID = ADT.SubjectID;

sensitivity(:,1) = scanAngleDegrees2./greenQuantaPerSecLambdaIntoEye;
sensitivity(:,2) = scanAngleDegrees2./redQuantaPerSecLambdaIntoEye;
tbl.sensitivity = sensitivity;

% minimize sum of squared residuals between V_lambda and data

fun = @(a) sum((a.*V_Lambda_quantal(wls==543) - tbl.sensitivity(:,1)).^2 + (a.*V_Lambda_quantal(wls==680) - tbl.sensitivity(:,2)).^2);
x0 =  median(tbl.sensitivity(:,1))./V_Lambda_quantal(wls==543);
x = fminsearch(fun, x0);

%plot(wls, log10(x.*V_Lambda_quantal), 'LineWidth', 1.5, 'Color', 'k');


% do same minimization as above for custom v lambda
S = [380 1 401];
wls = SToWls(S);

sensitivityWls = [543 680];
coneSpectralPeaks = [555.5 530.3 420.7; 558.9 530.3 420.7; 563.4 530.3 420.7];
offset = 5;
for subj = 1:numel(subjects)
    subjSensitivity = sensitivity(strcmpi(ADT.SubjectID, subjects{subj}));
    params = getSubjectParams(subjects{subj});
    for p = 1:size(coneSpectralPeaks,1)
        params.nomogram.lambdaMax = transpose(coneSpectralPeaks(p,:));
    %     unbleachedCustomPhotoreceptors = customConeSensitivities(params);
    %     [fractionBleachedFromIsom, ~, projectorRetIrradianceQuantaPerLambdaUm2Sec] = ...
    % projectorConeBleach(whiteRGB, projectorLum, projectorParams, eyeParams, unbleachedCustomPhotoreceptors);
    %     params.fractionPigmentBleached.value = fractionBleachedFromIsom;
    %     customPhotoreceptors{subj,p} = customConeSensitivities(params);

    photoreceptors{subj,p} = photoreceptorSensitivitiesFromScratch(params);
        
        L = photoreceptors{subj,p}.isomerizationAbsorptanceBleached(1,:);
        M = photoreceptors{subj,p}.isomerizationAbsorptanceBleached(2,:);
        V_lambda_unscaled = params.propL.*L + (1-params.propL).*M;
        V_lambda_fitted = fitVlambda(wls, V_lambda_unscaled, sensitivityWls, subjSensitivity, 543);
        plot(wls + offset.*(-1)^subj, log10(V_lambda_fitted), 'Color', subjColors(subj,:), 'LineStyle', lineStyles{p}, 'LineWidth', 1.5);


    end
end


% time for boxplots 

lambdaTest = [680 543];
boxWidth = 12;
lineWidth = 1.5;
offset = 5;
for wl = 1:2
    for subj = 1:numel(subjects)
        subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
        b{wl, subj} = boxchart(sensitivityWls(wl).*ones(size(sensitivity(subjIdx,wl))) + offset.*(-1)^subj, log10(sensitivity(subjIdx, wl)),...
            'BoxFaceColor', subjColors(subj,:), 'BoxEdgeColor', subjColors(subj,:), 'BoxWidth', boxWidth, 'LineWidth', lineWidth, 'MarkerColor', subjColors(subj,:));

    end
end

%set(gca, 'yscale', 'log')
xlim([513 710])

set(gcf, 'Units', 'centimeters')
set(gcf, 'Position', [0.5   0.5   9   9])

set(gca, 'Position', [0.2 0.2 0.75 0.75])

set(gca, 'fontsize', axFontSize, 'TickDir', tickDir, 'LineWidth', defaultLineWidth, 'TickLength', tickLength)
ylabel('Log Sensitivity (quanta^{-1} s deg^2)', 'fontsize', 14)
xlabel('Wavelength (nm)', 'FontSize', 14)
xticks([543 680]);

function V_lambda_fitted = fitVlambda(V_lambda_wls, V_lambda_unscaled, sensitivity_wls, sensitivity, fitWls)

[~,ia] = intersect(V_lambda_wls, fitWls); 
[~, ib] = intersect(sensitivity_wls, fitWls);

V_lambda_at_fitWls = V_lambda_unscaled(ia);

if numel(V_lambda_at_fitWls,1) ~= size(sensitivity,2)
    V_lambda_at_fitWls = transpose(V_lambda_at_fitWls);
else
end

fun = @(k) sum((k.*V_lambda_at_fitWls - sensitivity(:,ib)).^2, 'all');
x0 =  mean(median(sensitivity(:,ib),1)./V_lambda_at_fitWls');
x = fminsearch(fun, x0);

V_lambda_fitted = x.*V_lambda_unscaled;
