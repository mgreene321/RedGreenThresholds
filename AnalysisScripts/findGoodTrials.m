function goodDeliveries = findGoodTrials(dataFullFile, saveFlag, useDeliveryCentroids)
if nargin < 3
    useDeliveryCentroids = 0;
end

if nargin < 2
    saveFlag = 1;
end

if nargin < 1
    filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
    dataPath = uigetdir(filepath);
    [fileName, pathName] = uigetfile(dataPath);
    dataFullFile = fullfile(pathName, fileName);
    load(dataFullFile);
else
    [dataPath, dataFile, fext] = fileparts(dataFullFile);
    load(dataFullFile);
end

expDate = char(regexp(dataPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match')); % gets mm_dd_yyyy
expTime = char(regexp(dataPath , '\d{1,2}_\d{1,2}_\d{1,2}$', 'Match'));
%expTime = expTime{1};

subjectId = regexp(dataPath,  '\d{5}[LR]|\d{5}.[LR]', 'Match');
subjectId = erase(char(subjectId), filesep);

% Determine redGreenthresholds path

root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
saveFolder = fullfile(subjectFolder, [expDate], expTime);

if ~saveFlag
    if isfile(fullfile(saveFolder, 'goodDeliveries.mat')) && isfile(fullfile(saveFolder, 'findAllCrossesOutput.mat'))
        load(fullfile(saveFolder, 'goodDeliveries.mat'))
        load(fullfile(saveFolder, 'findAllCrossesOutput.mat'))
    else
        saveFlag = 1;
    end
end

if saveFlag
    expParams = exp_data.expParams;
    %stimDiam = expParams.stimDiam;
    %stimDur = expParams.stimDur;
    tca = [expParams.redTCA_X expParams.redTCA_Y expParams.greenTCA_X expParams.greenTCA_Y];
%     vidDurFrames = 30;
%     
%     % Search for cross within frames of interest, which are the frames during
%     % which the cross should nominally appear +/ 4 frames
%     stimFrames = (vidDurFrames/2 - floor(stimDur/2)) : (vidDurFrames/2 + floor(stimDur/2));
%     framesOfInterest = (stimFrames(1) - 4):(stimFrames(end) + 4);
    
    % Determine the frames of each video for which the cross is found, get all
    % frames for which cross is found, and average all frames for which cross
    % is found
    %[cross_XY,~] = findAllCrosses(framesOfInterest, dataPath, 0);
    
    findAllCrossesOutput = importdata(fullfile(saveFolder, 'findAllCrossesOutput.mat'));
    cross_XY = findAllCrossesOutput.cross_XY;
    % Analysis done in coordinate system of 712 x 712 stabilized frames
    
    if ~useDeliveryCentroids
        targetX = expParams.targetX + 100; targetY = expParams.targetY + 100;
    else
        delivered_X = cross_XY(:,1,1); delivered_Y = cross_XY(:,2,1);
        outliers_X = abs(delivered_X - mean(delivered_X, 'omitnan'))>= 1.5*std(delivered_X, 'omitnan');
        outliers_Y = abs(delivered_Y - mean(delivered_Y, 'omitnan'))>= 1.5*std(delivered_Y, 'omitnan');
        targetX = mean(delivered_X(~outliers_X), 'omitnan'); targetY = mean(delivered_Y(~outliers_Y), 'omitnan');
    end
    
    % Stimulus offsets
    [X Y] = meshgrid(0:expParams.gridHeight - 1, 0:expParams.gridWidth - 1);
    stimOff = expParams.stimDiam .* [X(:) Y(:)];
    numLocs = expParams.gridWidth * expParams.gridWidth;
    stimLocs(:,:,1) = stimOff + [targetX targetY]; %IR
    stimLocs(:,:,2) = stimLocs(:,:,1) + tca(1:2); % red
    stimLocs(:,:,3) = stimLocs(:,:,1) +  tca(3:4); %green
    
    data = exp_data.data_matrix;
    goodFrames = false(length(cross_XY), 1);
    
    % There may be some cases where the stimulus cross in one channel is hidden
    % by the cross in another channel. Let's identify and fix those cases
    % first.
    
    for i = 1:length(cross_XY)
        trialNum = cross_XY(i,4,1);
        loc = data(trialNum, 1);
        aom = data(trialNum, 3);
        stimLoc = stimLocs(loc,:, aom + 1);
        
        if size(data, 2) == 6
            stimDur = expParams.stimDur;
            stimDiam = expParams.stimDiam;
        else
            stimDiam = expParams.stimDiam(data(trialNum,4));
            stimDur = expParams.stimDur(data(trialNum,5));
            
        end
        
        
        if any(isnan(cross_XY(i,1:2,aom+1))) && ...
                sum(all(stimLocs(loc,:,:) - stimLoc == [0 0])) >= 2 % another cross occurs where cross of interest should be
            
            ChannelsWithOverlappingCrosses = find(all(stimLocs(loc,:,:) - stimLoc == [0 0]));
            ChannelsWithOverlappingCrosses = find(ChannelsWithOverlappingCrosses ~= aom + 1); %here 1 = IR, 2 = red, 3 = green
            if ~isempty(ChannelsWithOverlappingCrosses)
                channel = ChannelsWithOverlappingCrosses(1); % = aom + 1
                tco = stimLoc - stimLocs(loc,:,channel);
                cross_XY(i,1:2, aom+1) = cross_XY(i,1:2, channel) + tco;
            end
        end
    end
    
    delivery_tol_arcmin = 0.75;
    PPD = mean([exp_data.sessionParams.PPD_X exp_data.sessionParams.PPD_Y]);
    delivery_tol_pix = PPD * delivery_tol_arcmin / 60;
    
    % check whether delivery was good
    for i = 1:length(cross_XY)
        trialNum = cross_XY(i,4,1);
        
        loc = data(trialNum, 1);
        aom = data(trialNum, 3);
        
          if size(data, 2) == 6
            stimDur = expParams.stimDur;
            stimDiam = expParams.stimDiam;
        else
            stimDiam = expParams.stimDiam(data(trialNum,4));
            stimDur = expParams.stimDur(data(trialNum,5));
            
        end
        
        stimLoc = stimLocs(loc,:, aom + 1);

        % Check if cross is within ~0.5 arcmin of where we expect it to be
        if abs(cross_XY(i, 1:2, aom + 1) - stimLoc) <= round(delivery_tol_pix)
            goodFrames(i) = true;
        else
            goodFrames(i) = false;
        end
        
        %If IR cross isn't there, it's a bad frame
        
        if any(isnan(cross_XY(i,1:2,1)))
            goodFrames(i) = false;
        else
        end 
    end
    
    %cross_XY = cross_XY(goodFrames,:,:);
    %vidsWithCrosses =  unique(cross_XY(goodFrames,4,:)); % numbers of videos where crosses have been found
    
    % Check that stimulus was on for criterion number of frames
    
    cross_XY_good = cross_XY(goodFrames,:,:);
    
    %For each trial/video
    trialNums = unique(cross_XY_good(:,4,1));
    
    goodDeliveries =[];
    numGoodTrials = 0;
    for t = 1:numel(trialNums)
        
        trialNum = trialNums(t);
        aom = data(trialNum, 3);
        
        
          if size(data, 2) == 6
            stimDur = expParams.stimDur;
            stimDiam = expParams.stimDiam;
        else
            stimDiam = expParams.stimDiam(data(trialNum,4));
            stimDur = expParams.stimDur(data(trialNum,5));
            
        end
        
        
        crossData = cross_XY_good(cross_XY_good(:,4,1) == trialNum,:,:);
         ntupletIsThere = diff(crossData(:,3,aom+1), stimDur - 1) == 0;

        
        % check that there are n consecutive frames, n = stimDur
        if sum(ntupletIsThere) >= 1
            %just one good triplet of consecutive frames
            
            N = size(crossData,1);
            midPointIdx = round(median(1:N));
            temp = crossData(midPointIdx-(round(stimDur-1)/2):midPointIdx+(round(stimDur-1)/2),:,:);
            numGoodTrials = numGoodTrials + 1;
        else
            % multiple triplets, choose the first
            
            temp = [];
        end
        goodDeliveries = vertcat(goodDeliveries, temp);
        
    end
    
    if isempty(goodDeliveries) % possible that stimuli were persistently shifted to another location
        goodDeliveries = findGoodTrials(dataFullFile, 1,1);
        
    else
         %save average cross locations
    delivered_X = goodDeliveries(:,1,1); delivered_Y = goodDeliveries(:,2,1);
    mean_delivered_X = mean(delivered_X, 'omitnan'); mean_delivered_Y = mean(delivered_Y, 'omitnan');
    [X Y] = meshgrid(0:expParams.gridHeight - 1, 0:expParams.gridWidth - 1);
    stimOff = expParams.stimDiam .* [X(:) Y(:)];
    meanDeliveryLocs = stimOff + [mean_delivered_X mean_delivered_Y]; %IR
    
    if ~isdir(saveFolder)
        mkdir(saveFolder)
    end
    
    save(fullfile(saveFolder,'goodDeliveries.mat'), 'goodDeliveries')
    
    % Save info about delivery statistics
    numTrials = exp_data.expParams.numStair * ...
        exp_data.expParams.trialsPerStair * ...
        exp_data.expParams.gridHeight * ...
        exp_data.expParams.gridWidth * 2;
    
    proportionGoodTrials = numGoodTrials/max(trialNums);
    
    stdX = std(delivered_X, 'omitnan');
    stdY = std(delivered_Y, 'omitnan');
    
    stdX_arcmin = 60 * stdX / PPD;
    stdY_arcmin = 60 * stdY / PPD;
    
    save(fullfile(saveFolder, 'deliveryStatistics.mat'), 'proportionGoodTrials', 'numTrials', 'meanDeliveryLocs', 'stdX_arcmin', 'stdY_arcmin');
        
    end
    
else
end

end