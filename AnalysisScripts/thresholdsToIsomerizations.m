
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));

%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
pupilDiamMm = 7;
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

%% Set up photoreceptors structure

photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
photoreceptors.nomogram.lambdaMax = [563 530.3 420.7]';
photoreceptors = FillInPhotoreceptors(photoreceptors);

S = photoreceptors.nomogram.S;
wvl = SToWls(S);
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
photoreceptorsBleached.fractionPigmentBleached.value = projectorFractionBleachedFromIsom;
photoreceptorsBleached = FillInPhotoreceptors(photoreceptorsBleached);

%%
subjects = unique(analyzedDataTable.SubjectID);
% for each subject
for s = 1:length(subjects)
    locCounter = 1;
    subjectId = subjects{s};
    subjectTable = analyzedDataTable(strcmpi(analyzedDataTable.SubjectID,subjectId),:);
    folders = unique(subjectTable.Folder);
    % for each folder
    for f = 1:length(folders)
        folderTable = subjectTable(strcmpi(folders{f}, subjectTable.Folder),:);
        blocks = unique(folderTable.Block);
        % for each block (1 analysisOuput folder per block)
        for b = 1:length(blocks)
            block = blocks(b);
            % load analysisOutput folder and get red and green powers, thresholds in au
            analysisOutput = importdata(fullfile(RedGreenThresholdsPath,subjectId, folders{f}, ['analysisOutput' num2str(block) '.mat']));
            redPowerWattsAt680nm = analysisOutput.redPowerWattsAt680nm;
            greenPowerWattsAt543nm = analysisOutput.greenPowerWattsAt543nm;
            redThresholdPowerAU = analysisOutput.redThresholdPowerAU;
            greenThresholdPowerAU = analysisOutput.greenThresholdPowerAU;
            
            % for each location, estimate effective quanta
            
            for i = 1:length(analysisOutput.greenThresholdPowerAU)
                tRedAU = analysisOutput.redThresholdPowerAU(i);
                tGreenAU = analysisOutput.greenThresholdPowerAU(i);
                
                [~, tRedRU] = convertAOArbitraryToRealUnits(tRedAU,redPowerWattsAt680nm,680, 22,scanAngleDeg, photoreceptorsBleached);
                
                % incident red quanta
                incidentRedQuanta{s}(locCounter,:) = sum(tRedRU.rawCornIrradianceQuantaPerLambdaUm2SecIn) .* stimRenderTimeSec .* pupilAreaMm2 .* 1e6;
                transmittedRedQuanta{s}(locCounter,:) = sum(tRedRU.transmittedQuantaPerLambdaUm2Sec) .* stimRenderTimeSec .* pupilAreaMm2 .* 1e6;
                effectiveRedQuantaLMS{s}(locCounter,:) = sum(tRedRU.effectiveQuantaPerLambdaUm2Sec,2) * stimRenderTimeSec * ISareaUm2;
                
                [~, tGreenRU] = convertAOArbitraryToRealUnits(tGreenAU,greenPowerWattsAt543nm,543, 22,scanAngleDeg, photoreceptorsBleached);
                incidentGreenQuanta{s}(locCounter,:) = sum(tGreenRU.rawCornIrradianceQuantaPerLambdaUm2SecIn) * stimRenderTimeSec  * pupilAreaMm2 * 1e6;%* stimAreaUm2;
                transmittedRedQuanta{s}(locCounter,:) = sum(tGreenRU.transmittedQuantaPerLambdaUm2Sec) .* stimRenderTimeSec .* pupilAreaMm2 .* 1e6;
                effectiveGreenQuantaLMS{s}(locCounter,:) = sum(tGreenRU.effectiveQuantaPerLambdaUm2Sec,2) * stimRenderTimeSec * ISareaUm2;
                
                locCounter = locCounter + 1;
                
            end
        end
    end
end