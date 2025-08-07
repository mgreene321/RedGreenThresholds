% Validate cross finding
warning('on');
filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
[vidName, vidPath] = uigetfile(fullfile(filepath, '*stabilized.avi'));

vidNum = str2double(regexp(vidName, '(?<=_)\d{3}', 'match'));


expDate = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match')); % gets mm_dd_yyyy
expTime = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{1,2}(?=\\)', 'Match'));

fprintf([expDate '\\' expTime '\\' num2str(vidNum) '\n']);
subjectId = char(regexp(vidPath,  '\d{5}[LR]|\d{5}.[LR]', 'Match'));
subjectId = erase(subjectId, filesep);

root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);

load(fullfile(subjectFolder, expDate, expTime, 'findAllCrossesOutput.mat'))
vidObj = VideoReader(fullfile(vidPath, vidName));
vidFrames = squeeze(read(vidObj));

crosses = cross_XY(cross_XY(:,4,1) == vidNum,:,:);
frameNums = crosses(:,3,1);

h = figure; hold on
imshow(vidFrames(:,:,1)), axis equal, xlim([1 size(vidFrames,2)]);
set(gca, 'YDir', 'reverse');
title(['Frame ' num2str(frameNums(1))])

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
        if f - 1 < min(frameNums)
            warning('Cannot rewind more');
        else
            f = f - 1;
            imshow(vidFrames(:,:, f)); hold on
            
            IR_cross = crosses(crosses(:,3,1) == f,1:2,1);
            red_cross = crosses(crosses(:,3,2) == f,1:2,2);
            green_cross = crosses(crosses(:,3,3) == f,1:2,3);
            plot(IR_cross(1), IR_cross(2), 'ko', 'MarkerFaceColor', 'w')
            plot(red_cross(1), red_cross(2), 'ro', 'MarkerFaceColor', 'r')
            plot(green_cross(1), green_cross(2), 'go', 'MarkerFaceColor', 'g')
            
            
            title(['Frame ' num2str(f)]);
        end
    elseif strcmpi(btn, 'rightarrow')
        if f + 1 > max(frameNums)
            warning('Cannot fast forward more');
        else
            f = f + 1;
            imshow(vidFrames(:,:, f)); hold on
             IR_cross = crosses(crosses(:,3,1) == f,1:2,1);
            red_cross = crosses(crosses(:,3,2) == f,1:2,2);
            green_cross = crosses(crosses(:,3,3) == f,1:2,3);
            plot(IR_cross(1), IR_cross(2), 'ko', 'MarkerFaceColor', 'w')
            plot(red_cross(1), red_cross(2), 'ro', 'MarkerFaceColor', 'r')
            plot(green_cross(1), green_cross(2), 'go', 'MarkerFaceColor', 'g')
            
            title(['Frame ' num2str(f)]);
        end
    elseif strcmpi(btn, 'return') || strcmpi(btn, 'escape')
        browseVideo = 0;
    else
    end
    
    
    
end