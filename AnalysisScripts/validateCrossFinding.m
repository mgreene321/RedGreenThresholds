% Validate cross finding


%% Load in data table
vidPath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
RawDataTable = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));

%% Select video
warning('on');
filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
[vidName, vidPath] = uigetfile(fullfile(filepath, '*stabilized.avi'), 'Select a video');

vidNum = str2double(regexp(vidName, '(?<=_)\d{3}', 'match'));
expDate = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match')); % gets mm_dd_yyyy
expTime = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{1,2}(?=\\)', 'Match'));

fprintf([expDate '\\' expTime '\\' num2str(vidNum) '\n']);


%% Determine subject, experiment date folder
subjectId = char(regexp(vidPath,  '\d{5}[LR]|\d{5}.[LR]', 'Match'));
subjectId = erase(subjectId, filesep);
idx = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.ExperimentDate, expDate) & strcmpi(RawDataTable.ExperimentTime, expTime);
expDateFolder = unique(RawDataTable.Folder(idx));


%% Load experiment data
dataDir = dir(fullfile(vidPath, 'RedGreenThresholds*.mat'));
temp = importdata(fullfile(dataDir(1).folder, dataDir(1).name));
exp_data = temp.exp_data;
tca = [exp_data.expParams.redTCA_X...
    exp_data.expParams.redTCA_Y...
    exp_data.expParams.greenTCA_X...
    exp_data.expParams.greenTCA_Y];
%% Determine row for this trial in RawDataTable

rowIdx = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.Folder, expDateFolder{1})...
    & strcmpi(RawDataTable.ExperimentTime, expTime) & RawDataTable.Trial == vidNum;
%% Determine stimulus channel

aom = RawDataTable.Channel(rowIdx);
if aom == 1
    
    tco = tca(1:2);
elseif aom == 2
    
    tco = tca(3:4);
end

stimColors = {'r', 'g'};
%% Load cross locations file

subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
load(fullfile(subjectFolder, expDateFolder{1}, expTime, 'findAllCrossesOutput.mat'))

%% Load in delivery analysis
load(fullfile(subjectFolder, expDateFolder{1}, expTime, 'goodDeliveries.mat'));
goodTrials = unique(goodDeliveries(:,4,1));

if ismember(vidNum, goodTrials)
    titleFragment = ' Good delivery';
else
    titleFragment = ' Bad delivery';
end



%% Load in daily image

meanFrame = importdata(fullfile(subjectFolder, expDateFolder{1}, 'TransferFromOCT', 'meanFrame.mat'));

%% Load in OCT image and cone location data

coneData = importdata(fullfile(subjectFolder, 'OCTMaster', 'coneData.mat'));
dailyData = importdata(fullfile(subjectFolder, expDateFolder{1}, 'TransferFromOCT', 'dailyData.mat'));

OCTImage = coneData.Raw_image_fundus_scaled;

%% Load cone apertures
coneApertures = importdata(fullfile(subjectFolder, expDateFolder{1}, 'TransferFromOCT', 'coneApertures.mat'));
coneApertures = coneApertures./max(coneApertures(:));
coneApertures = coneApertures(:,:,1:3);
%% Figure out daily and master correspondence points

dailyConeLocs = str2double(dailyData(:, [1 2]));
badIdx = any(isnan(dailyConeLocs),2) | any(dailyConeLocs < 0,2);
masterConeLocs = str2double(coneData.coneData(:, [1 2]));
masterCorrPts = masterConeLocs(~badIdx,:);
dailyCorrPts = dailyConeLocs(~badIdx,:);


%% Read in video
vidObj = VideoReader(fullfile(vidPath, vidName));
vidFrames = squeeze(read(vidObj));
crosses = cross_XY(cross_XY(:,4,1) == vidNum,:,:);
frameNums = crosses(:,3,1);

%% Figure out cross locations in master image

[xm, ym] = transformImageCoordinates(reshape(crosses(:,1,:), [], 1, 1), reshape(crosses(:,2,:), [],1,1), 1, [], dailyCorrPts, masterCorrPts);
%masterCrosses =reshape([xm ym], [],2,3);
masterCrosses(:,1,:) = reshape(xm, [], 1, 3);
masterCrosses(:,2,:) =reshape(ym, [], 1, 3);


%% Allow user to scroll through video
h = figure; hold on

% Plot video frames in first column
subplot(2,2,1); hold on
title('Stabilized video')
imagesc(vidFrames(:,:,1)),colormap gray, axis equal, axis off, xlim([1 size(vidFrames,2)]), ylim([1 size(vidFrames,1)]);
set(gca, 'YDir', 'reverse');
sgtitle(['Frame ' num2str(frameNums(1)) titleFragment])

% Plot daily image
subplot(2,2,2); hold on
title('Mean image')
imagesc(meanFrame), colormap gray, axis equal, axis off,  xlim([1 size(meanFrame,2)]),ylim([1 size(meanFrame,1)]) ;
set(gca, 'YDir', 'reverse');

% Plot OCT master image
subplot(2,2,3); hold on
title('OCT Master Image')
imagesc(OCTImage), colormap gray, axis equal, axis off,  xlim([1 size(OCTImage, 2)]), ylim([1 size(OCTImage, 1)]);
set(gca, 'YDir', 'reverse');

% Plot cone apertures
subplot(2,2,4); hold on
title('Cone spectral map')
imagesc(coneApertures), axis equal, axis off, xlim([1 size(coneApertures,2)]), ylim([1 size(coneApertures,1)]);
set(gca, 'YDir', 'reverse');
browseVideo = 1;
KbName('UnifyKeyNames');
set(h, 'KeyPressFcn', 'uiresume')

f = frameNums(1);
while browseVideo
    uiwait;
    
    try
        btn = get(h, 'CurrentKey');
    catch
        browseVideo = 0;
    end
    
    % disp(btn);
    if strcmpi(btn, 'leftarrow')
        
        try
            delete(rectHand);
        catch
        end
        if f - 1 < min(frameNums)
            warning('Cannot rewind more');
        else
            f = f - 1;
            sgtitle(['Frame ' num2str(f) titleFragment]);
            IR_cross = crosses(crosses(:,3,1) == f,1:2,1);
            red_cross = crosses(crosses(:,3,2) == f,1:2,2);
            green_cross = crosses(crosses(:,3,3) == f,1:2,3);
            
            % Mark crosses with white, red, or green marker according to
            % channel
            subplot(2,2,1), hold on; colormap gray, axis equal, axis off, xlim([1 size(vidFrames,2)]), ylim([1 size(vidFrames,1)]);
            imagesc(vidFrames(:,:, f)); hold on
            
            
            
            plot(IR_cross(1), IR_cross(2), 'ko', 'MarkerFaceColor', 'w')
            plot(red_cross(1), red_cross(2), 'ro', 'MarkerFaceColor', 'r')
            plot(green_cross(1), green_cross(2), 'go', 'MarkerFaceColor', 'g')
            
            
            % Plot marker for channel that was on, on meanFrame
            subplot(2,2,2); imagesc(meanFrame), colormap gray, axis equal, axis off,  xlim([1 size(meanFrame,2)]),ylim([1 size(meanFrame,1)]) ;
            stimX =crosses(crosses(:,3,1) == f, 1, aom+1); stimY = crosses(crosses(:,3,1) == f, 2, aom+1);
            plot(stimX, stimY, 'Marker', 'o', 'MarkerEdgeColor', stimColors{aom}, 'MarkerFaceColor', stimColors{aom});
            
            % PLot on master
            
            subplot(2,2,3); imagesc(OCTImage), colormap gray, axis equal, axis off,  xlim([1 size(OCTImage, 2)]), ylim([1 size(OCTImage, 1)]);
            masterStimX = masterCrosses(crosses(:,3,1) ==f, 1, aom+1); masterStimY = masterCrosses(crosses(:,3,1) == f, 2, aom+1);
            plot(masterStimX, masterStimY, 'Marker', 'o', 'MarkerEdgeColor', stimColors{aom}, 'MarkerFaceColor', stimColors{aom});
            % Plot cone apertures
            subplot(2,2,4); hold on
            
            
            set(gca, 'YDir', 'reverse');
            % Draw stim
            try
                imagesc(coneApertures), axis equal, axis off, xlim([1 size(coneApertures,2)]), ylim([1 size(coneApertures,1)]);
                rectHand = drawrectangle('Position', [stimX-tco(1)-11, stimY-tco(2)-11, 21, 21], 'LineWidth', 1.5, 'Color', stimColors{aom}, 'FaceAlpha', 0.3, 'InteractionsAllowed', 'none');
                title(['L: ' num2str(RawDataTable.L(rowIdx)) ' M: ' num2str(RawDataTable.M(rowIdx)) 'S: ' num2str(RawDataTable.S(rowIdx))]);
                plot(stimX, stimY, 'yo');
            catch
                imagesc(coneApertures), axis equal, axis off, xlim([1 size(coneApertures,2)]), ylim([1 size(coneApertures,1)]);
            end
            
            
        end
    elseif strcmpi(btn, 'rightarrow')
        try
            delete(rectHand);
        catch
        end
        if f + 1 > max(frameNums)
            warning('Cannot fast forward more');
        else
            f = f + 1;
            sgtitle(['Frame ' num2str(f) titleFragment]);
            IR_cross = crosses(crosses(:,3,1) == f,1:2,1);
            red_cross = crosses(crosses(:,3,2) == f,1:2,2);
            green_cross = crosses(crosses(:,3,3) == f,1:2,3);
            
            subplot(2,2,1),imagesc(meanFrame), colormap gray, axis equal, axis off,  xlim([1 size(meanFrame,2)]),ylim([1 size(meanFrame,1)]) ;
            imagesc(vidFrames(:,:, f)); hold on; colormap gray, axis equal, axis off, xlim([1 size(vidFrames,2)]), ylim([1 size(vidFrames,1)]);
            
            plot(IR_cross(1), IR_cross(2), 'ko', 'MarkerFaceColor', 'w')
            plot(red_cross(1), red_cross(2), 'ro', 'MarkerFaceColor', 'r')
            plot(green_cross(1), green_cross(2), 'go', 'MarkerFaceColor', 'g')
            
            
            
            % Plot marker for channel that was on, on meanFrame
            subplot(2,2,2)
            imagesc(meanFrame), colormap gray, axis equal, axis off,  xlim([1 size(meanFrame,2)]),ylim([1 size(meanFrame,1)]) ;
            stimX =crosses(crosses(:,3,1) == f, 1, aom+1); stimY = crosses(crosses(:,3,1) == f, 2, aom+1);
            plot(stimX, stimY, 'Marker', 'o', 'MarkerEdgeColor', stimColors{aom}, 'MarkerFaceColor', stimColors{aom});
            
            
            % PLot on master
            
            subplot(2,2,3), imagesc(OCTImage), colormap gray, axis equal, axis off,  xlim([1 size(OCTImage, 2)]), ylim([1 size(OCTImage, 1)]);
            masterStimX = masterCrosses(crosses(:,3,1) ==f, 1, aom+1); masterStimY = masterCrosses(crosses(:,3,1) == f, 2, aom+1);
            plot(masterStimX, masterStimY, 'Marker', 'o', 'MarkerEdgeColor', stimColors{aom}, 'MarkerFaceColor', stimColors{aom});
            
            % Plot cone apertures
            subplot(2,2,4); hold on
            
            
            set(gca, 'YDir', 'reverse');
            % Draw stim
            try
                imagesc(coneApertures), axis equal, axis off, xlim([1 size(coneApertures,2)]), ylim([1 size(coneApertures,1)]);
                rectHand = drawrectangle('Position', [stimX-tco(1)-11, stimY-tco(2)-11, 21, 21], 'LineWidth', 1.5, 'Color', stimColors{aom}, 'FaceAlpha', 0.3, 'InteractionsAllowed', 'none');
                title(['L: ' num2str(RawDataTable.L(rowIdx)) ' M: ' num2str(RawDataTable.M(rowIdx)) 'S: ' num2str(RawDataTable.S(rowIdx))]);
                plot(stimX, stimY, 'yo');
                
            catch
                imagesc(coneApertures), axis equal, axis off, xlim([1 size(coneApertures,2)]), ylim([1 size(coneApertures,1)]);
            end
            
        end
    elseif strcmpi(btn, 'return') || strcmpi(btn, 'escape')
        browseVideo = 0;
    else
    end
end
