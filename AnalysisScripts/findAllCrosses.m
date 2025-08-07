function [cross_XY, meanFrame] = findAllCrosses(dataPath, saveFlag)
% INPUTS
% - framesOfInterest: indices of video frames in which to search for cross

% OUTPUTS
% - allCrosses (uint8 array): each row corresponds to a frame. First row
%   element is cross x-coordinate, second row element is cross y-coordinate,
%   third row element is frame index, fourth row element is video index.
%
%- allFrames (uint8 array): First dimension is frame row, second
%   dimension si frame column, third dimension is frame index, fourht
%   dimension is video index.
%
% - meanFrame: average all allFrames


vidDurFrames = 30; %hardcoded: number of frames per video

% Save outputs (overwriting any old ones) by default
if nargin < 2
    saveFlag = 1;
end

% If folder containing videos to analyze isn't specified, let user select
% folder
if nargin < 1
    filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
    dataPath = uigetdir(filepath); %full file path leading to where videos are saved
end


% Load the experiment data mat file (why wasn't I doing this before???)
RGDataFileDir = dir(fullfile(dataPath, 'RedGreenThresholds*.mat'));
RGDataFileName = RGDataFileDir.name;
RedGreenThresholdsData = importdata(fullfile(dataPath,RGDataFileName));
exp_data = RedGreenThresholdsData.exp_data;
data_matrix = exp_data.data_matrix;
stimDurs = exp_data.expParams.stimDur;
stimDiams = exp_data.expParams.stimDiam;


% Extract the experiment date, time, and subject from the video folder name
expDate = regexp(dataPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match'); % gets mm_dd_yyyy
expDate = expDate{1};

expTime = regexp(dataPath , '\d{1,2}_\d{1,2}_\d{1,2}$', 'Match');
expTime = expTime{1};

subjectId = regexp(dataPath,  '\d{5}[LR]|\d{5}.[LR]', 'Match');
subjectId = erase(subjectId{1}, filesep);

% Find the project folder, which is where everything gets saved
root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
saveFolder = fullfile(subjectFolder, expDate, expTime);

% If we aren't saving, load saved outputs if they exist. Otherwise,
% generate new outputs and save them
if ~saveFlag
    if isfile(fullfile(saveFolder, 'findAllCrossesOutput.mat'))
        load(fullfile(saveFolder, 'findAllCrossesOutput.mat'))
    else
        saveFlag = 1;
    end
end

if saveFlag
    % Extract video numbers and determine number of videos
    vidDir = dir(fullfile(dataPath, '*_stabilized.avi'));
    vidNames = {vidDir.name};
    vidNumArray = cellfun(@(x) regexp(x, '\d*(?=_stabilized.avi)', 'match'), vidNames);
    vidNums = cellfun(@(x) str2double(x), vidNumArray);
    numVids = length(vidNums);
    
    % Make stimulus cross
    cross_diam = 11;
    cross = zeros(cross_diam);
    cross(6,:) = 1;
    cross(:,6) = 1;
    
    % for each video, and for each frame, find cross
    
    i = 1; % frame counter
    % framesToIgnore = false(length(frame_vids),1);
    
    sumFrames = 0;
    pixelCumSum = false(712,712);
    
    for vidIdx = 1:numVids
        
        try
            vidObj = VideoReader(fullfile(dataPath,vidDir(vidIdx).name));
            if size(data_matrix, 2) == 6 % old format data matrix where size and duration columns arent included
                stimDur = 3;
                stimDiam = 21;
            else
                stimDiam = stimDiams(data_matrix(vidNums(vidIdx), 4));
                stimDur = stimDurs(data_matrix(vidNums(vidIdx), 5));
            end
            
            % Search for cross within frames of interest, which are the frames during
            % which the cross should nominally appear +/- 4 frames
            stimFrames = (vidDurFrames/2 - floor(stimDur/2)) : (vidDurFrames/2 + floor(stimDur/2));
            framesOfInterest = (stimFrames(1) - 4):(stimFrames(end) + 4);
            framesOfInterest(framesOfInterest < 1 | framesOfInterest > vidDurFrames) = [];
            
            for frameIdx = 1:numel(framesOfInterest)
                
                frameNum = framesOfInterest(frameIdx);
                frame = read(vidObj, frameNum);
                
                temp_IR_frame = frame .* uint8(frame == 255);
                temp_red_frame = frame .* uint8(frame == 253);
                temp_green_frame = frame .* uint8(frame == 251);
                
                filtered_IR_frame = filter2(cross, temp_IR_frame);
                filtered_red_frame = filter2(cross, temp_red_frame);
                filtered_green_frame = filter2(cross, temp_green_frame);
                
                [y_IR, x_IR] = find(filtered_IR_frame == max(filtered_IR_frame(:)));
                [y_red, x_red] = find(filtered_red_frame == max(filtered_red_frame(:)));
                [y_green, x_green] = find(filtered_green_frame == max(filtered_green_frame(:)));
                
                % For each channel, confirm that a cross is present at the peak of
                % the filtered frame. If it is, set row of the appropriate layer (layer 1 = IR, layer 2  = red, layer 3 = green) of
                % the matrix cross_XY equal to the cross x and y coordinates, and
                % the frame and video in which the cross was found.
                
                % Generate black image
                blank_im = zeros([size(frame) 3]);
                
                if crossMightBeThere(x_IR, y_IR)
                    x_IR = x_IR(1); y_IR = y_IR(1);
                    if x_IR-floor(cross_diam/2) > 0 && y_IR-floor(cross_diam/2) > 0
                        % Draw IR cross
                        blank_im(y_IR, x_IR-floor(cross_diam/2):x_IR+floor(cross_diam/2),1) = 1;
                        blank_im(y_IR-floor(cross_diam/2):y_IR+floor(cross_diam/2), x_IR,1) = 1;
                    else
                        x_IR = nan; y_IR = nan;
                    end
                else
                    x_IR = nan; y_IR =nan;
                end
                
                if crossMightBeThere(x_red, y_red)
                    x_red = x_red(1); y_red = y_red(1);                   
                    if x_red-floor(cross_diam/2) > 0 && y_red-floor(cross_diam/2) > 0
                        % Draw red cross
                        
                        blank_im(y_red, x_red-floor(cross_diam/2):x_red+floor(cross_diam/2),2) = 1;
                        blank_im(y_red-floor(cross_diam/2):y_red+floor(cross_diam/2), x_red,2) = 1;
                        
                    else
                        x_red = nan; y_red = nan;
                    end
                else
                    x_red = nan; y_red = nan;
                end
                
                
                if crossMightBeThere(x_green, y_green)
                    x_green = x_green(1); y_green = y_green(1);
           
                    
                    % Draw green cross
                    if x_green-floor(cross_diam/2) > 0 && y_green-floor(cross_diam/2) > 0
                        blank_im(y_green, x_green-floor(cross_diam/2):x_green+floor(cross_diam/2),3) = 1;
                        blank_im(y_green-floor(cross_diam/2):y_green+floor(cross_diam/2), x_green,3) = 1;
                    else
                        x_green = nan; y_green = nan;
                    end
                else
                    x_green = nan; y_green = nan;
                end
                
                % Handle overlapping crosses
                
                
                hypoIRCross = xor(blank_im(:,:,1), sum(blank_im,3)>1) .* blank_im(:,:,1);
                hypoRedCross = xor(blank_im(:,:,2), sum(blank_im,3)>1) .*blank_im(:,:,2);
                hypoGreenCross = xor(blank_im(:,:,3), sum(blank_im,3)>1) .* blank_im(:,:,3);
                
                % Make hypothetical crosses thicker
                avgFilter = ones(3,3);
                IRCrossSearch = conv2(hypoIRCross, avgFilter, 'same') > 0;
                RedCrossSearch = conv2(hypoRedCross, avgFilter, 'same') > 0;
                GreenCrossSearch = conv2(hypoGreenCross, avgFilter, 'same') > 0;
                
                if sum(hypoIRCross(:)) > 1 && crossIsThere(x_IR, y_IR, frame, hypoIRCross, IRCrossSearch, 255)
                    cross_XY(i,:,1) = [x_IR y_IR frameNum vidNums(vidIdx)];
                else
                    cross_XY(i,:,1) = [nan nan frameNum vidNums(vidIdx)];
                end
                
                
                %RED
                
                if sum(hypoRedCross(:)) > 1 && crossIsThere(x_red, y_red, frame, hypoRedCross, RedCrossSearch, 253)
                    cross_XY(i,:,2) = [x_red y_red frameNum vidNums(vidIdx)];
                else
                    cross_XY(i,:,2) = [nan nan frameNum vidNums(vidIdx)];
                end
                
                
                %GREEN
                
                if sum(hypoGreenCross(:)) > 1 && crossIsThere(x_green, y_green, frame, hypoGreenCross, GreenCrossSearch, 251)
                    cross_XY(i,:,3) = [x_green y_green framesOfInterest(frameIdx) vidNums(vidIdx)];
                else
                    cross_XY(i,:,3) = [nan nan frameNum vidNums(vidIdx)];
                end
                
                
                if all(all(isnan(cross_XY(i,1:2,:))))
                    framesToIgnore(i) = true;
                else
                    
                    sumFrames = sumFrames + double(frame);
                    pixelCumSum = pixelCumSum + (frame > 0);
                end
                
                i = i + 1;
                
            end
        catch
        end
    end
    
    %meanFrame = sumFrames./sum(~framesToIgnore);
    meanFrame = sumFrames./pixelCumSum;
    
    
    if ~isdir(saveFolder)
        mkdir(saveFolder)
    end
    
    save(fullfile(saveFolder,'findAllCrossesOutput.mat'), 'cross_XY', 'meanFrame')
    
end

% Auxiliary functions

function bool = crossMightBeThere(x_cross, y_cross)

% When filtering image for cross, two nearby peaks may be found

if numel(x_cross) == 1 && numel(y_cross) == 1
    bool = true;
elseif numel(x_cross) == 2 && numel(y_cross) == 2
    
    if abs(diff(x_cross)) <= 1 && abs(diff(y_cross)) <= 1
        bool = true;
    else
        bool = false;
    end
else
    bool = false;
    
end

function bool = crossIsThere(x_cross, y_cross, frame, hypoCross, crossSearch, crossVal)

% Confirms the presence of a cross if either the entire horizontal bar is
% found, the entire vertical bar is found, or if the two halves of the
% veritcal bar are found slightly offset from one another.

maskedFrame = uint8(crossSearch).*frame;

% conditions to find cross

searchRegionFilled = sum(maskedFrame(:) == crossVal) >= sum(hypoCross(:)) - 1;
horizontalBarFound = sum(maskedFrame(y_cross,:) == crossVal) >= sum(hypoCross(y_cross,:)) - 1;
verticalBarFound =  sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal, 'all') >= sum(hypoCross(:,x_cross)) - 1 &&...
    ~any(sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal, 2) > 1); %entire, possibly sheared vertical bar found


horizontalBarNotTooTall = sum(sum(maskedFrame(y_cross-1:y_cross+1, :) == crossVal,1) > 2) <= 4;
verticalBarNotTooWide = sum(sum(maskedFrame(:, x_cross-1:x_cross+1) == crossVal,2) > 2) <= 4;

if horizontalBarFound && horizontalBarNotTooTall
    bool = true;
elseif verticalBarFound && verticalBarNotTooWide
    bool = true;
else
    bool = false;
end
