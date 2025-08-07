close all;
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
%Load analyzed data table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);
% load raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);

propL = RDT.L./(RDT.L + RDT.M);

inflectionProp = 0.5;
hetero = 1-(inflectionProp^-1).*abs(propL-inflectionProp);
% distanceFromGlobalPropL = 1-(globalPropL.^-1).*abs(propL-globalPropL);


indices = getLocationIndicesInRawDataTable;

for i = 1:size(RDT,1)
    temp = cellfun(@(x) find(x), indices, 'UniformOutput', false);
    ADT_idx(i) = find(cellfun(@(x) ismember(i, x), temp));

    if RDT.Channel(i) == 1
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.RedThresholdPowerAU(ADT_idx(i));

    elseif RDT.Channel(i) == 2
        propThresh(i,:) = RDT.IntensityAU(i)./ADT.GreenThresholdPowerAU(ADT_idx(i));
    end
end

% s cone proximity
[~, meanDist, SconesUnderStim] = SconeProximity(1, 1, 0);


expectedResponse = double(RDT.Channel == RDT.Color);
% tbl = table(RDT.SubjectID, ADT_idx', propThresh, hetero, RDT.Channel, meanDist', expectedResponse);
% tbl.Properties.VariableNames = {'SubjectID', 'Location', 'Intensity', 'Heterogeneity', 'Channel', 'SconeDistance', 'ExpectedResponse'};

subjectIndices = cellfun(@(x) find(strcmpi(x,subjects)), RDT.SubjectID);

tbl = table(subjectIndices, ADT_idx', propThresh, hetero, RDT.Channel, meanDist', expectedResponse);
tbl.Properties.VariableNames = {'SubjectIndex', 'Location', 'Intensity', 'Heterogeneity', 'Channel', 'S', 'ExpectedResponse'};

tbl = tbl(RDT.YesNo == 1,:);

% center variables
intensityRef = log2(1);
heteroRef = 0; %2/3;
subjRef = min(tbl.SubjectIndex);
channelRef = min(tbl.Channel);
SRef = mean(tbl.S);

tbl.Intensity = log2(tbl.Intensity);
tbl.Heterogeneity = tbl.Heterogeneity - heteroRef;
tbl.SubjectIndex = tbl.SubjectIndex - subjRef;
tbl.Channel = tbl.Channel - channelRef;
tbl.S = tbl.S - SRef;


% as in paper
% modelspec_full_hetero= 'ExpectedResponse ~ Intensity*Heterogeneity + (1|SubjectIndex) + (1|Channel)';
% modelspec_full_scone= 'ExpectedResponse ~ Intensity*SconeDistance + (1|SubjectIndex) + (1|Channel)';
%
%
%
% %modelspec_reduced = 'ExpectedResponse ~ Intensity + (1|SubjectID) + (1|Channel)';
% %%%%%
%
% %modelspec_full = 'ExpectedResponse ~ Heterogeneity + Intensity:Heterogeneity + Intensity*SconeDistance + (1|SubjectID) + (1|Channel)';
% %modelspec_full = 'ExpectedResponse ~ SubjectIndex + Channel  + Intensity*Heterogeneity + Intensity:SconeDistance';
% modelspec_reduced = 'ExpectedResponse ~ SubjectIndex + Channel  + Intensity*Heterogeneity';
%
% glme_full = fitglme(tbl, modelspec_full, 'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');

mdl_glm_full = 'ExpectedResponse ~ SubjectIndex+ Channel + Intensity + Heterogeneity + S';
mdl_glme_full = 'ExpectedResponse ~ Intensity + Heterogeneity + S + (1|SubjectIndex) + (1|Channel)';
mdl_glme_reduced =  'ExpectedResponse ~ Intensity + S + (1|SubjectIndex) + (1|Channel)';

mdl_FEglme = 'ExpectedResponse ~ Intensity + Heterogeneity + S';

glm_full = fitglm(tbl, mdl_glm_full, Distribution="binomial");
glme_full = fitglme(tbl, mdl_glme_full,'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');
FEglme = fitglme(tbl, mdl_FEglme, 'Distribution','Binomial', 'Link', 'Logit', 'FitMethod', 'Laplace');
glme_reduced = fitglme(tbl, mdl_glme_reduced,'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');

%glm_reduced = fitglm(tbl, modelspec_reduced, Distribution="binomial");
%glme_reduced = fitglme(tbl, modelspec_reduced, 'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');
% glmeGlobal = fitglme(tblGlobal, modelspecGlobal, 'Distribution', 'Binomial', 'Link', 'Logit');

% glme_hetero = fitglme(tbl, modelspec_full_hetero, 'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');
% glme_scone = fitglme(tbl, modelspec_full_scone, 'Distribution', 'Binomial',  'Link', 'Logit', 'FitMethod', 'Laplace');

% Compare model predictions with empirical data
y = predict(glm_full, tbl);

heteroEdges = [0 quantile(tbl.Heterogeneity + heteroRef, 2) 1];
heteroEdges = round(2*heteroEdges,1)/2; % round to nearest 0.5
SEdges = [min(tbl.S+SRef) quantile(tbl.S + SRef,2) max(tbl.S+SRef)];
intEdges = -2:5;


%% P vs. hetero
[N_hetero, ~, binHetero] = histcounts(tbl.Heterogeneity + heteroRef, heteroEdges);
[N_hetero_exp, ~, binHetero_exp] = histcounts(tbl.Heterogeneity(tbl.ExpectedResponse == 1) + heteroRef, heteroEdges);

% median hetero value in each bin
uniqueBinHetero = unique(binHetero);
for i = 1:numel(uniqueBinHetero)
    medianHetero(i) = median(tbl.Heterogeneity(binHetero == uniqueBinHetero(i)) + heteroRef);
end
P_hetero = N_hetero_exp./N_hetero;

figure; hold on
plot(tbl.Heterogeneity + heteroRef, y, 'bo');
plot(medianHetero(uniqueBinHetero>0), P_hetero, 'ro-', 'MarkerFaceColor','r');

%% P vs. intensity

[N_int, ~, binInt] = histcounts(tbl.Intensity, intEdges);
[N_int_exp, ~, binInt_exp] = histcounts(tbl.Intensity(tbl.ExpectedResponse == 1), intEdges);

% median hetero value in each bin
uniqueBinInt = unique(binInt);
for i = 1:numel(uniqueBinInt)
    medianInt(i) = median(tbl.Intensity(binInt == uniqueBinInt(i)));
end
P_int = N_int_exp./N_int;
figure; hold on
plot(tbl.Intensity, y, 'bo');
plot(medianInt(uniqueBinInt>0), P_int, 'ro-', 'MarkerFaceColor','r');

%% P vs. S

[N_S, ~, binS] = histcounts(tbl.S + SRef, SEdges);
[N_S_exp, ~, binS_exp] = histcounts(tbl.S(tbl.ExpectedResponse == 1), SEdges);

% median hetero value in each bin
uniqueBinS = unique(binS);
for i = 1:numel(uniqueBinS)
    medianS(i) = median(tbl.S(binS == uniqueBinS(i)) + SRef);
end
P_S = N_S_exp./N_S;
figure; hold on
plot(tbl.S + SRef, y, 'bo');
plot(medianS(uniqueBinS>0), P_S, 'ro-', 'MarkerFaceColor','r');


%% P vs. subject

figure; hold on;
uniqueSubjectIndex = unique(tbl.SubjectIndex);
for i = 1:numel(uniqueSubjectIndex)
    P_subj_obs(i) = sum(tbl.ExpectedResponse(tbl.SubjectIndex == uniqueSubjectIndex(i)))./numel(tbl.ExpectedResponse(tbl.SubjectIndex==uniqueSubjectIndex(i)));
    P_subj_pred(i) = mean(y(tbl.SubjectIndex == uniqueSubjectIndex(i)));
end
temp = transpose([P_subj_obs; P_subj_pred]);
bar(categorical(subjects), temp);


%% P vs. channel

figure; hold on;
uniqueChannel = unique(tbl.Channel);
for i = 1:numel(uniqueChannel)
    P_channel_obs(i) = sum(tbl.ExpectedResponse(tbl.Channel == uniqueChannel(i)))./numel(tbl.ExpectedResponse(tbl.Channel==uniqueChannel(i)));
    P_channel_pred(i) = mean(y(tbl.Channel == uniqueChannel(i)));
end
temp = transpose([P_channel_obs; P_channel_pred]);
bar(categorical({'680 nm', '543 nm'}), temp);

%% P vs. intensity, grouped by hetero
[N, ~, ~, binInt, binHetero] = histcounts2(tbl.Intensity, tbl.Heterogeneity+heteroRef, intEdges, heteroEdges);
[N_exp, ~, ~, binInt_exp, binHetero_exp] = histcounts2(tbl.Intensity(tbl.ExpectedResponse == 1), tbl.Heterogeneity(tbl.ExpectedResponse == 1)+heteroRef, intEdges, heteroEdges);
P = N_exp./N;

uniqueBinInt = unique(binInt);
uniqueBinHetero = unique(binHetero);


for i = 1:numel(uniqueBinInt)
    medianIntByHetero(i) = median(tbl.Intensity(binInt == uniqueBinInt(i)));
end

figure; hold on
for i = 1:size(P,2)
    plot(medianIntByHetero(uniqueBinInt>0), P(:,i), 'ro-')
    plot(tbl.Intensity(binHetero == uniqueBinHetero(i)), y(binHetero == uniqueBinHetero(i)), 'o')
end



%% P vs intensity, grouped by hetero for each channel and subject
colors = [1 0 0; 0 1 0];
markers = {'^', 's', 'o'};
for subj = 1:numel(uniqueSubjectIndex) % for each subject
    figure; hold on; title(['Grouped by heterogeneity, ' subjects{subj}]);
     xlabel('log2_ threshold multiple');
    ylabel('P(chromatic response)');
      xlim([-2 5]); ylim([0 1]);
      set(gca, 'FontSize', 12)
    for ch = 1:numel(uniqueChannel) % for each channel
      
        tempTbl = tbl(tbl.SubjectIndex == uniqueSubjectIndex(subj) & tbl.Channel == uniqueChannel(ch),:);
        tempy = y(tbl.SubjectIndex == uniqueSubjectIndex(subj) & tbl.Channel == uniqueChannel(ch));

        [N, ~, ~, binInt, binHetero] = histcounts2(tempTbl.Intensity, tempTbl.Heterogeneity+heteroRef, intEdges, heteroEdges);
        [N_exp, ~, ~, binInt_exp, binHetero_exp] = histcounts2(tempTbl.Intensity(tempTbl.ExpectedResponse == 1),...
            tempTbl.Heterogeneity(tempTbl.ExpectedResponse == 1)+heteroRef, intEdges, heteroEdges);

        N = N.*(N > 10);
        P = N_exp./N;

        uniqueBinInt = unique(binInt);
        uniqueBinHetero = unique(binHetero);

        medianIntByHetero = [];
        for i = 1:numel(uniqueBinInt)
            medianIntByHetero(i) = median(tempTbl.Intensity(binInt == uniqueBinInt(i)));
        end

        for i = 1:numel(uniqueBinHetero)

            if uniqueBinHetero(i) > 0
            plot(tempTbl.Intensity(binHetero == uniqueBinHetero(i)), tempy(binHetero == uniqueBinHetero(i)), 'o')
            else
            end
        end
        for i = 1:size(P,2)
            
            plot(medianIntByHetero(uniqueBinInt>0), P(uniqueBinInt(uniqueBinInt>0),i),'Color', 'k',...
                'Marker', markers{i}, 'MarkerSize', 9, 'MarkerFaceColor', colors(ch,:), 'MarkerEdgeColor', 'k', 'LineStyle', '-', 'LineWidth', 2)
            
        end


    end
end

%% P vs intensity, grouped by S for each channel and subject
colors = [1 0 0; 0 1 0];
for subj = 1:numel(uniqueSubjectIndex) % for each subject
    figure; hold on; title(['Grouped by S cone distance, ' subjects{subj}]);
    xlim([-2 5]); ylim([0 1]);
    xlabel('log_2 threshold multiple');
    ylabel('P(chromatic response)')
    set(gca, 'fontsize', 12)
    for ch = 1:numel(uniqueChannel) % for each channel
      
        tempTbl = tbl(tbl.SubjectIndex == uniqueSubjectIndex(subj) & tbl.Channel == uniqueChannel(ch),:);
        tempy = y(tbl.SubjectIndex == uniqueSubjectIndex(subj) & tbl.Channel == uniqueChannel(ch));

        [N, ~, ~, binInt, binS] = histcounts2(tempTbl.Intensity, tempTbl.S + SRef, intEdges, SEdges);
        [N_exp, ~, ~, binInt_exp, binS_exp] = histcounts2(tempTbl.Intensity(tempTbl.ExpectedResponse == 1),...
            tempTbl.S(tempTbl.ExpectedResponse == 1)+SRef, intEdges, SEdges);

        N = N.*(N > 10);
        P = N_exp./N;

        uniqueBinInt = unique(binInt);
        uniqueBinS = unique(binS);

        medianIntByS = [];
        for i = 1:numel(uniqueBinInt)
            medianIntByS(i) = median(tempTbl.Intensity(binInt == uniqueBinInt(i)));
        end

        for i = 1:numel(uniqueBinS)

            if uniqueBinS(i) > 0
            plot(tempTbl.Intensity(binS == uniqueBinS(i)), tempy(binS == uniqueBinS(i)), 'o')
            else
            end
        end
        for i = 1:size(P,2)
            
            plot(medianIntByS(uniqueBinInt>0), P(uniqueBinInt(uniqueBinInt>0),i),'Color', 'k',...
                'Marker', markers{i}, 'MarkerSize', 9, 'MarkerFaceColor', colors(ch,:),'MarkerEdgeColor', 'k' ,'LineStyle', '-', 'LineWidth', 2)
            
        end


    end
end

