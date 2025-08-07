function [thresholdsPowerAU, psiParamsFit, goodData] = RedGreenThresholds_Analysis(redPower, greenPower)
addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\mQUESTPlus-master'))
%convert from uW to W
redPower = redPower/1e6;
greenPower = greenPower/1e6;
% Select files for analysis

tempFolder = {'Z:\Local_Share\Imaging_Data\AOVIS'};
selectFiles = 1;
i = 1;
while selectFiles
    [tempFile, tempFolder] = uigetfile(fullfile(tempFolder, '*.mat'), 'Select data file');
    selectFiles = ischar(tempFile);
    if ~selectFiles; break; end
    fileNames{i,:} = tempFile;
    dataPaths{i,:} = tempFolder;
    i = i + 1;
end

expDate = regexp(dataPaths , '\d{1,2}_\d{1,2}_\d{4}', 'Match'); % gets mm_dd_yyyy
expTime = regexp(dataPaths, '\d{1,2}_\d{1,2}_\d{1,2}(?=\\)', 'Match');
subjectId = regexp(dataPaths,  '\d{5}[LR]|\d{5}.[LR]', 'Match');

expDate = [expDate{:}]';
expTime = [expTime{:}]';

subjectId = cellfun(@(x) erase(x, filesep), subjectId, 'UniformOutput', false); % get rid of \
subjectId = [subjectId{:}]';

if length(unique(expDate)) > 1
    warning('Data from multiple days selected')
end

if length(unique(subjectId)) > 1
    warning('Data from multiple subjects selected')
end


% naming output file: check if there are output files so far

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId{1});
saveFolder = fullfile(subjectFolder, [expDate{1}]);

analysisDir = dir(fullfile(saveFolder, 'analysisOutput*'));

analysisFileName = fullfile(saveFolder,['analysisOutput' num2str(size(analysisDir,1)+1) '.mat']);


% Compile data, cross locations, frames

for i = 1:length(fileNames)
    folderName = dataPaths{strcmp(fileNames, fileNames{i,:}),:};
    
    %goodDeliveries = findGoodTrials(fullfile(folderName, fileNames{i,:}),0);
    goodDeliveries = importdata(fullfile(saveFolder, expTime{i}, 'goodDeliveries.mat'));
    tempGoodTrialNums = unique(goodDeliveries(:,4,1));
    %tempGoodTrialNums = findGoodTrials(fullfile(folderName, fileNames{i,:}),0);
    
    goodTrialNums{i,:} = tempGoodTrialNums;
    load(fullfile(folderName, fileNames{i,:}));
    data = exp_data.data_matrix;
    goodData{i,:} = data(tempGoodTrialNums,:);
    
    
end

goodData = vertcat(goodData{:});
psiParamsFit = fitAndPlotPFs(goodData);
thresholdsPowerAU = computeThresholds(psiParamsFit, goodData, redPower, greenPower);


% analysisOutput.redThresholdPowerAU = thresholdsPowerAU(:,1);
% analysisOutput.greenThresholdPowerAU = thresholdsPowerAU(:,2);

analysisOutput.redThresholdPowerAU = cellfun(@(x) x(:,1), thresholdsPowerAU, 'UniformOutput', false);
analysisOutput.greenThresholdPowerAU = cellfun(@(x) x(:,2), thresholdsPowerAU, 'UniformOutput', false);


% added 4/9/24
analysisOutput.pixelsPerSecond = 20e6;
analysisOutput.redPowerWattsAt680nm = redPower;
analysisOutput.greenPowerWattsAt543nm = greenPower;
anaysisOutput.pixelsPerFrame = 512*512;
analysisOutput.stimDurFrames = exp_data.expParams.stimDur;
analysisOutput.stimDiamPix = exp_data.expParams.stimDiam;

analysisOutput.psiParamsFit = psiParamsFit;
analysisOutput.goodData = goodData;
analysisOutput.expTime = expTime;

save(analysisFileName, 'analysisOutput')

function psiParamsFit = fitAndPlotPFs(goodData)

locIdx = unique(goodData(:,1));
colorIdx = unique(goodData(:,3));
numLocs = length(locIdx);

if size(goodData,2) == 6
    numDiams = 1;
    numDurs = 1;
else
    numDurs = length(unique(goodData(:,5)));
    numDiams = length(unique(goodData(:,4)));
end

%figure; hold on;
lineColors = {'r', 'g'};
lineColors = 0.*[0.5 0.5 0.5; 0 0 0];
markerColors =[1 0 0; 0 1 0];
lineStyles = {'-', ':'};
markerShapes = {'o', 's'};

%fit weibulls
for diam = 1:numDiams
    for dur = 1:numDurs
        figure; hold on
        for l = 1:numLocs
            loc = locIdx(l);
            subplot(ceil(sqrt(numLocs)),floor(sqrt(numLocs)),loc); hold on
            for c = 1:2
                color = colorIdx(c);
                
                if size(goodData,2) == 6
                tempData = goodData(goodData(:,1) == loc & goodData(:,3) == color,:);
                else
                     tempData = goodData(goodData(:,1) == loc & goodData(:,3) == color & goodData(:,4) == diam & goodData(:,5) == dur,:);
                end
                
                tempData(tempData(:,end-2) == 0,:) =[];
                
                trialData = [];
                for i = 1:length(tempData)
                    trialData(i).stim = 20.*log10(tempData(i,end-2));
                    trialData(i).outcome = tempData(i,end-1) + 1;
                end
                % Set up fake quest+
                qpPF = @qpPFWeibull;
                nOutcomes = 2;
                stimParamsDomainList = {20.*log10(0:0.001:1)}; %for cum norm
                psiParamsDomainList  = {20.*log10(0:0.001:1), 3, 0.01, 0.01}; %mean, SD are fixed
                stopRule = 'nTrials';
                q = qpInitialize('qpPF', qpPF,...
                    'nOutcomes', nOutcomes,...
                    'stimParamsDomainList', stimParamsDomainList,...
                    'psiParamsDomainList', psiParamsDomainList,...
                    'stopRule', stopRule, 'noentropy', true);
                for i = 1:length(tempData)
                    q = qpUpdate(q, trialData(i).stim, trialData(i).outcome);
                end
                psiParamsIndex = qpListMaxArg(q.posterior);
                psiParamsQuest = q.psiParamsDomain(psiParamsIndex,:);
                psiParamsFit{l,c, diam, dur} = qpFit(q.trialData, q.qpPF, psiParamsQuest, q.nOutcomes,...
                    'lowerBounds', [20.*log10(0) 0.5 0.01 0.01], 'upperBounds', [20.*log10(1) 20 0.01 0.01]);
                
                % struct with two fields: stim (stimulus intensity) and
                % outcomeCounts (first column gives number of "not seen" responses,
                % second column gives number of "seen" responses, for correspodning
                % stimulus intensity)
                stimCounts = qpCounts(qpData(q.trialData), q.nOutcomes);
                stim = [stimCounts.stim];
                stimFine = 20.*log10(linspace(0,1,1000))';
                plotProportionsFit = qpPFWeibull(stimFine, psiParamsFit{l,c});
                
                  plot(10.^(stimFine./20),plotProportionsFit(:,2),'LineStyle', lineStyles{c},'Color',lineColors(c,:),'LineWidth',3);
                % for each stimulus intensity
                for cc = 1:length(stimCounts)
                    nTrials(cc) = sum(stimCounts(cc).outcomeCounts);
                    pCorrect(cc) = stimCounts(cc).outcomeCounts(2)/nTrials(cc);
                end
%                 for cc = 1:length(stimCounts)
%                     h = scatter(10.^(stim(cc)./20),pCorrect(cc),100,'Marker', markerShapes{c},'MarkerEdgeColor','none','MarkerFaceColor',lineColors{c},...
%                         'MarkerFaceAlpha',nTrials(cc)/max(nTrials),'MarkerEdgeAlpha',nTrials(cc)/max(nTrials));
%                 end
              
%                 %             xlabel('Stimulus Value');
%                 %             ylabel('Proportion Correct');
                
                % plot it with binning that makes more sense
                
                temp = [stimCounts.outcomeCounts];
                temp = transpose(reshape(temp, 2,[]));
                %[N, edges, bins] = histcounts(stim, 'NumBins', 20, 'BinMethod', 'fd');
                [N, edges, bins] = histcounts(stim, 'NumBins', 20);
                uniqueBins = unique(bins);
                for i = 1:numel(unique(bins))
                    bin = uniqueBins(i);
                    pCorr(i) = sum(temp(bins == bin,2))/sum(nTrials(bins == bin));
                end
                for i = 1:numel(unique(bins))
                    bin = uniqueBins(i);
                   
                    % normalize to max(N), the number of data in most
                    % populous bin
                    markerSize = 20.*sqrt(sum(bins == uniqueBins(i))/max(N));
                   % plot(10.^(mean(stim(bins==uniqueBins(i)))/20), pCorr(i), 'Marker', 'o', 'MarkerFaceColor', lineColors{c}, 'MarkerEdgeColor', [0.5 0.5 0.5],'LineWidth', 2, 'MarkerSize', markerSize)
                    plot(10.^(mean(stim(bins==uniqueBins(i)))/20), pCorr(i), 'Marker',...
                        'o', 'MarkerFaceColor', markerColors(c,:), 'MarkerEdgeColor', lineColors(c,:),...
                        'LineWidth', 2, 'MarkerSize', markerSize, 'lineStyle', 'none')
                end
                
                
            end
        end
    end
end

% get 50% threshold

function thresholdsPowerAU= computeThresholds(psiParamsFit,goodData, redPower, greenPower)
locIdx = unique(goodData(:,1));
colorIdx = unique(goodData(:,3));
numLocs = length(locIdx);

if size(goodData,2) == 6
    numDiams = 1;
    numDurs = 1;
else
    numDurs = length(unique(goodData(:,5)));
    numDiams = length(unique(goodData(:,4)));
end

FoS = 0.5;
stimFine = 20.*log10(linspace(0,1,1000))';
for l = 1:numLocs
    for c = 1:2
        for diam = 1:numDiams
            for dur = 1:numDurs
                try
                    params = psiParamsFit{l,c, diam, dur};
                    
                    psiFun = qpPFWeibull(stimFine, psiParamsFit{l,c, diam, dur});
                    psiFun = psiFun(:,2);
                    
                    % restrict psiFun to domain over which it is changing
                    
                    unique_psiFun = unique(psiFun);
                    count_psiFun = hist(psiFun, unique_psiFun);
                    singleton_psiFun = unique_psiFun(count_psiFun == 1);
                    singleton_stimFine = stimFine(count_psiFun == 1);
                    psiFunInv = @(x) interp1(singleton_psiFun, singleton_stimFine, x);
                    
                    thresholdsPowerAU{diam, dur}(l,c) = 10.^(psiFunInv(FoS)./20);
                catch
                end
            end
        end
    end
end

% pix_per_sec = 20e6; pix_per_frame = 512*512; numFrames = 3;
% thresholds_Energy_au = thresholds_Power_au .* (1/pix_per_sec) .* pix_per_frame .* numFrames;
%
% redThresholds_Energy_J = redPower .* thresholds_Energy_au(:,1);
% greenThresholds_Energy_J = greenPower .* thresholds_Energy_au(:,2);
%
% thresholds = [redThresholds_Energy_J greenThresholds_Energy_J];

