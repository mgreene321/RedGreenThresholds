function [sproximity, meanDist, SconesUnderStim] = SconeProximity(byTrialFlag, nearestCones, outsideStim)
PPD = 560;
% S cone proximity analysis
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
%Load analyzed data table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
% Load in raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);

previousConeDataFile = '';
if byTrialFlag
    for i = 1:size(RDT,1)
        subjectId = RDT.SubjectID{i};
        expFolder = RDT.Folder{i};

        stimX = RDT.Xc(i);
        stimY = RDT.Yc(i);

        coneDataFile = fullfile(RedGreenThresholdsPath, subjectId, expFolder,'TransferFromOCT', 'dailyData.mat');
        missingConesFile = fullfile(RedGreenThresholdsPath, subjectId, expFolder,'TransferFromOCT', 'likelyMissingCones.mat');

        if ~strcmpi(coneDataFile, previousConeDataFile)
            coneData = importdata(coneDataFile);
            likelyMissingCones = importdata(missingConesFile);
            cone_locs = str2double(coneData(:, 1:2));
            cone_labels = coneData(:,3);

            badIdx = any(cone_locs < 0,2) | any(isnan(cone_locs),2);

            cone_locs = cone_locs(~badIdx,:);
            cone_labels = cone_labels(~badIdx);

            all_cone_locs = [cone_locs; likelyMissingCones];
            unique_cone_locs = unique(all_cone_locs, 'rows');
            cone_spacing = estimateConeSpacing(unique_cone_locs, PPD, 'triangulation');

        else
        end

        previousConeDataFile = coneDataFile;

        Scone_locs = cone_locs(cone_labels=="S",:);

        % figure out cone locs within stimulated area
        xv = stimX + RDT.StimDiamPix(i).*[-1/2 1/2 1/2 -1/2 -1/2];
        yv = stimY + RDT.StimDiamPix(i).*[-1/2 -1/2 1/2 1/2 -1/2];
        in = inpolygon(Scone_locs(:,1), Scone_locs(:,2), xv, yv);

        if outsideStim
        D = pdist2([stimX stimY], Scone_locs(~in,:));
        else
            D = pdist2([stimX stimY], Scone_locs);
        end

        D = D.*60/PPD;
        D = D./cone_spacing;

        sproximity(i) = sum(1./(1+D));

        D = sort(D);
        meanDist(i) = mean(D(nearestCones));

        SconesUnderStim(i) = sum(in);


    end

else
    for i = 1:size(ADT,1)
        subjectId = ADT.SubjectID{i};
        expFolder = ADT.Folder{i};

        stimX = ADT.meanXc(i);
        stimY = ADT.meanYc(i);

        coneDataFile = fullfile(RedGreenThresholdsPath, subjectId, expFolder,'TransferFromOCT', 'dailyData.mat');
        missingConesFile = fullfile(RedGreenThresholdsPath, subjectId, expFolder,'TransferFromOCT', 'likelyMissingCones.mat');

        if ~strcmpi(coneDataFile, previousConeDataFile)
            coneData = importdata(coneDataFile);
            likelyMissingCones = importdata(missingConesFile);
            cone_locs = str2double(coneData(:, 1:2));
            cone_labels = coneData(:,3);

            badIdx = any(cone_locs < 0,2) | any(isnan(cone_locs),2);

            cone_locs = cone_locs(~badIdx,:);
            cone_labels = cone_labels(~badIdx);

            all_cone_locs = [cone_locs; likelyMissingCones];
            unique_cone_locs = unique(all_cone_locs, 'rows');
            cone_spacing = estimateConeSpacing(unique_cone_locs, PPD, 'triangulation');

        else
        end

        previousConeDataFile = coneDataFile;

        Scone_locs = cone_locs(cone_labels=="S",:);

         % figure out cone locs within stimulated area
        xv = stimX + ADT.StimDiamPix(i).*[-1/2 1/2 1/2 -1/2 -1/2];
        yv = stimY + ADT.StimDiamPix(i).*[-1/2 -1/2 1/2 1/2 -1/2];
        in = inpolygon(Scone_locs(:,1), Scone_locs(:,2), xv, yv);


         if outsideStim
        D = pdist2([stimX stimY], Scone_locs(~in,:));
        else
            D = pdist2([stimX stimY], Scone_locs);
        end

        D = D.*60/PPD;
        D = D./cone_spacing;
        SconeDensity = sum(cone_labels=="S")./numel(cone_labels);
        sproximity(i) = sum(1./(1+D))./SconeDensity;

        D = sort(D);
        meanDist(i) = mean(D(nearestCones));
        SconesUnderStim(i) = sum(in);
        % as a check that things are being done right, compute L and M cones
        % stimSize = 21;
        % stimVertices = [stimX stimY] + fix(stimSize/2)*[-1 -1; 1 -1; -1 1; 1 1];
        % [inStim, onStim] = inpolygon(cone_locs(:,1), cone_locs(:,2), stimVertices(:,1), stimVertices(:,2));
        % L_by_loc(i) = sum((inStim | onStim) & cone_labels=="L");
        % M_by_loc(i) = sum((inStim | onStim) & cone_labels=="M");
        % S_by_loc(i) = sum((inStim | onStim) & cone_labels=="S");



    end
end
end

% redPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMRed.mat');
% greenPrimaryTbl = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\PR650MeasurementsOfAOPrimaries\meanAOMGreen.mat');
%
% % was SplineCMF
% Red_energy = transpose(SplineSpd(redPrimaryTbl.wavelength, redPrimaryTbl.power, wvl));
% Green_energy = transpose(SplineSpd(greenPrimaryTbl.wavelength, greenPrimaryTbl.power, wvl));
%
% S = [380 1 401];
% wvl = SToWls(S);
% load T_xyz1931;
% T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
% clear T_xyz1931 S_xyz1931
%
%
% indices = getLocationIndicesInRawDataTable;
% propThresh = nan(size(RDT,1),1);
% for i = 1:size(RDT,1)
%
%     % find correspodning data point in ADT
%     temp = cellfun(@(x) find(x), indices, 'UniformOutput', false);
%     ADT_idx = cellfun(@(x) ismember(i, x), temp);
%
%     if RDT.Channel(i) == 1
%         lambdaMax = 680;
%         maxLaserPowerWatts = 0.4*unique(ADT.RedPowerWattsAt680nm(ADT_idx));
%         energySpectrum = Red_energy;
%         propThresh(i,:) = RDT.IntensityAU(i)./ADT.RedThresholdPowerAU(ADT_idx);
%
%     elseif RDT.Channel(i) == 2
%         lambdaMax = 543;
%         maxLaserPowerWatts = 0.4*unique(ADT.GreenPowerWattsAt543nm(ADT_idx));
%         energySpectrum = Green_energy;
%         propThresh(i,:) = RDT.IntensityAU(i)./ADT.GreenThresholdPowerAU(ADT_idx);
%     end
%
% end