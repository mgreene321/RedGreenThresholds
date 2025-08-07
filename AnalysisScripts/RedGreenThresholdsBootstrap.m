function RedGreenThresholdsBootstrap

addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\mQUESTPlus-master'))
%convert from uW to W
%redPower = redPower/1e6;
%greenPower = greenPower/1e6;
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
    
    goodDeliveries = findGoodTrials(fullfile(folderName, fileNames{i,:}),0);
    
    tempGoodTrialNums = unique(goodDeliveries(:,4,1));
    %tempGoodTrialNums = findGoodTrials(fullfile(folderName, fileNames{i,:}),0);
    
    goodTrialNums{i,:} = tempGoodTrialNums;
    load(fullfile(folderName, fileNames{i,:}));
    data = exp_data.data_matrix;
    goodData{i,:} = data(tempGoodTrialNums,:);
end

goodData = vertcat(goodData{:});
[paramsValues, thresh, LL, exitFlag, boostrap] = fitAndPlotPFs(goodData);

function [psiParamsFit, thresh, LL, exitFlag, boostrap] =  fitAndPlotPFs(goodData)
PF = @PAL_CumulativeNormal;
parametric = 0;
guess = 0.01;
lapse = 0.01;
paramsFree = [1 1 0 1];
PF = @PAL_CumulativeNormal;
niter = 1000; % for bootstrapping
thresholdPercent = 50;

locIdx = unique(goodData(:,1));
colorIdx = unique(goodData(:,3));
numLocs = length(locIdx);

figure; hold on;
colors = [1 0 0; 0 1 0];
markerShapes = {'o', 's'};

%fit weibulls

for l = 1:numLocs
    loc = locIdx(l);
    subplot(ceil(sqrt(numLocs)),floor(sqrt(numLocs)),loc); hold on
    for c = 1:2
        color = colorIdx(c);
        tempData = goodData(goodData(:,1) == loc & goodData(:,3) == color,:);
        
        stimLevels = tempData(:,4);
        numSeen = tempData(:,5);
        outOfNum = ones(size(tempData(:,5)));
        
        
        [paramsValues, thresh, LL, exitFlag, boostrap] = fitPSF_RG(stimLevels, numSeen, outOfNum, PF, thresholdPercent, paramsFree, guess, lapse, niter, parametric);
        psiParamsFit{l,c} = paramsValues;
        
        plist = linspace(min(stimLevels),max(stimLevels), 100);
        [~, ~, idx] = histcounts(stimLevels, plist);
        x = linspace(0,1,1e3);
        y = PF(psiParamsFit{l,c}, x);
        plot(x,y,  'Color', colors(c,:), 'LineWidth', 2);
        for h = 1:length(plist)
            pred(h) = mean(tempData(idx == h,5));
            if length(tempData(idx == h,5)) == 0
                msize = 1/3;
            else
                msize = (2/3) * numel(tempData(idx == h,5));
                plot(plist(h),pred(h),'o','Color', colors(c,:), 'LineWidth', 2, 'MarkerFaceColor',[0.75 0.75 0.75],'markersize',msize,'HandleVisibility','off'); hold on
            end
        end
       
        
        
        
    end
end





