function F = makeStimulusDemoVideo

% good video:

%'1_25_2024\14_19_4\16';

%% Load in data table
vidPath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
RawDataTable = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));


%% Select video
warning('on');
filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
[vidName_unstab, vidPath] = uigetfile(fullfile(filepath, '*.avi'), 'Select a video');

vidNum = str2double(regexp(vidName_unstab, '(?<=_)\d{3}', 'match'));
expDate = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match')); % gets mm_dd_yyyy
expTime = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{1,2}(?=\\)', 'Match'));

fprintf([expDate '\\' expTime '\\' num2str(vidNum) '\n']);

%% Determine subject, experiment date folder
subjectId = char(regexp(vidPath,  '\d{5}[LR]|\d{5}.[LR]', 'Match'));
subjectId = erase(subjectId, filesep);
idx = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.ExperimentDate, expDate) & strcmpi(RawDataTable.ExperimentTime, expTime);
expDateFolder = unique(RawDataTable.Folder(idx));
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);

%% Load in delivery analysis
goodDeliveries = importdata(fullfile(subjectFolder, expDateFolder{1}, expTime, 'goodDeliveries.mat'));
goodTrials = unique(goodDeliveries(:,4,1));

% get frame numbers

if ismember(vidNum, goodTrials)
else
    [vidName_unstab, vidPath] = uigetfile(fullfile(filepath, '*.avi'), 'Select a video');

    vidNum = str2double(regexp(vidName_unstab, '(?<=_)\d{3}', 'match'));
    expDate = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{4}', 'Match')); % gets mm_dd_yyyy
    expTime = char(regexp(vidPath , '\d{1,2}_\d{1,2}_\d{1,2}(?=\\)', 'Match'));

    fprintf([expDate '\\' expTime '\\' num2str(vidNum) '\n']);
end

vidName_stab = [vidName_unstab(1:end-4) '_stabilized.avi'];

frameNums = goodDeliveries(goodDeliveries(:,4,1) == vidNum, 3,1);
frameNums = frameNums -1;
%% Determine IR cross locations in video

IR_cross_XY_unstab = findIRCrossInSingleRawVideo(fullfile(vidPath, vidName_unstab), frameNums);
IR_cross_XY_unstab = IR_cross_XY_unstab - 1;

meanX_unstab = mean(IR_cross_XY_unstab(:,1));
meanY_unstab = mean(IR_cross_XY_unstab(:,2));



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
    lambda = 680;
elseif aom == 2
    lambda = 543;
end

stim_unstab = blurStim(21, 512, lambda);
stim_stab = blurStim(21, 712, lambda);

%% Determine response
resp = RawDataTable.Color(rowIdx);
propL = RawDataTable.L(rowIdx)./(RawDataTable.L(rowIdx) + RawDataTable.M(rowIdx));
responses = {'Red', 'Green', 'Achromatic', 'Not Seen'};

%% stabilized cross location

IR_cross_XY_stab = goodDeliveries(goodDeliveries(:,4,1) == vidNum,1:2,1);


meanX_stab = mean(IR_cross_XY_stab(:,1));
meanY_stab = mean(IR_cross_XY_stab(:,2));


%% Get cone data
%cone_data = impordata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\20217R\2_27_2024\TransferFromOCT\dailyData.mat')
cone_data = importdata(fullfile(RedGreenThresholdsPath, subjectId, expDate,  'TransferFromOCT', 'dailyData.mat'));
cone_locs = str2double(cone_data(:, 1:2));
cone_labels = cone_data(:,3);

coneApertures = importdata(fullfile(RedGreenThresholdsPath, subjectId, expDate,  'TransferFromOCT', 'coneApertures.mat'));
coneApertures = coneApertures(:,:,1:3)./max(coneApertures(:,:,1:3), [], 'all');
coneApertures = imgaussfilt(coneApertures+0.1667, 2.5);
%% test plot

vidObj_unstab = VideoReader(fullfile(vidPath,vidName_unstab));

vidObj_stab = VideoReader(fullfile(vidPath, vidName_stab)); 

% h = figure; hold on
% set(gca, 'ydir', 'reverse')

xc = 512/2; yc = 512/2;

meanStimLoc = mean(IR_cross_XY_unstab,1);
x_mean = meanStimLoc(1);
y_mean = meanStimLoc(2);

%while isvalid(h)
unstabilized_gif_fileName = fullfile(RedGreenThresholdsPath, subjectId, ['unstabilized_' expDate '_' expTime '_' num2str(vidNum) '.gif']);
stabilized_gif_fileName = fullfile(RedGreenThresholdsPath, subjectId, ['stabilized_' expDate '_' expTime '_' num2str(vidNum) '.gif']);
coneApertures_gif_fileName = fullfile(RedGreenThresholdsPath, subjectId, ['coneApertures_' expDate '_' expTime '_' num2str(vidNum) '.gif']);

dt = 1/6;
patchSize = 140;


frames_unstab = im2double(read(vidObj_unstab));
frames_unstab = repmat(frames_unstab, [1 1 3 1]);

frames_stab = im2double(read(vidObj_stab));
frames_stab = repmat(frames_stab, [1 1 3 1]);

frames_coneApertures = repmat(coneApertures, [1 1 1 vidObj_stab.NumFrames]);

for f = 1:numel(frameNums)
    stimLoc_unstab = [IR_cross_XY_unstab(f,1) IR_cross_XY_unstab(f,2)];
    stim_shifted_unstab = circshift(stim_unstab, [stimLoc_unstab(2) stimLoc_unstab(1)] - [xc yc]);
    frames_unstab(:,:, aom, frameNums(f)) = frames_unstab(:,:, aom, frameNums(f)) + stim_shifted_unstab;

    stimLoc_stab = [IR_cross_XY_stab(f,1) IR_cross_XY_stab(f,2)];
    stim_shifted_stab = circshift(stim_stab, [stimLoc_stab(2) stimLoc_stab(1)] - [712 712]/2);
    frames_stab(:,:, aom, frameNums(f)+1) = frames_stab(:,:, aom, frameNums(f)+1) + stim_shifted_stab;

    frames_coneApertures(:,:, aom, frameNums(f)+1) = frames_coneApertures(:,:,aom, frameNums(f)+1) + stim_shifted_stab;

end

for i = 1:vidObj_unstab.NumFrames

    frame_unstab = frames_unstab(:,:,:,i);
    frame_unstab_patch = frame_unstab(meanY_unstab - patchSize/2:meanY_unstab + patchSize/2,...
        meanX_unstab - patchSize/2:meanX_unstab + patchSize/2,:);

    frame_stab = frames_stab(:,:,:,i);
    frame_stab_patch = frame_stab(meanY_stab - patchSize/2:meanY_stab + patchSize/2,...
        meanX_stab - patchSize/2:meanX_stab + patchSize/2,:);

    frame_coneApertures = frames_coneApertures(:,:,:,i);
    frame_coneApertures_patch =  frame_coneApertures(meanY_stab - patchSize/2:meanY_stab + patchSize/2,...
        meanX_stab - patchSize/2:meanX_stab + patchSize/2,:);

    [A1, map1] = rgb2ind(frame_unstab_patch, 256);
    [A2, map2] = rgb2ind(frame_stab_patch, 256);
    [A3, map3] = rgb2ind(frame_coneApertures_patch, 256);

    if i < 8
        % imwrite(A1, map1, retina_map_gif_fileName, "gif", LoopCount = Inf, DelayTime = dt);
    elseif i == 8
        imwrite(A1, map1, unstabilized_gif_fileName, "gif", LoopCount = Inf, DelayTime = dt);
        imwrite(A2, map2, stabilized_gif_fileName, "gif", LoopCount = Inf, DelayTime = dt);
        imwrite(A3, map3, coneApertures_gif_fileName, "gif", LoopCount = Inf, DelayTime = dt);

    elseif i < 22
        imwrite(A1, map1,unstabilized_gif_fileName, "gif", WriteMode = "append",  DelayTime = dt);
        imwrite(A2, map2,stabilized_gif_fileName, "gif", WriteMode = "append",  DelayTime = dt);
        imwrite(A3, map3,coneApertures_gif_fileName, "gif", WriteMode = "append",  DelayTime = dt);

    end
end


stim = blurStim(21, 712, lambda);

X_stab = IR_cross_XY_stab(end,1);
Y_stab = IR_cross_XY_stab(end,2);

stimLoc_unstab = [meanX_stab meanY_stab];
    stim_shifted_unstab = circshift(stim, round([stimLoc_unstab(2) stimLoc_unstab(1)]) - [712 712]/2);

coneApertures_stim = coneApertures.*(0.25+stim_shifted_unstab);
coneApertures_patch = coneApertures(meanY_stab - patchSize/2:meanY_stab + patchSize/2,...
        meanX_stab - patchSize/2:meanX_stab + patchSize/2,:);



function stim = blurStim(stimDiamPix, imSize, lambda)

pupilDiamMm = 6.5;
ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05;

stimTemplate = zeros(imSize, imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;

PSF = generate_PSF(ap_field, psf_pupil, zernike_pupil, field_size, lambda/1000, diff_limited, defocus);
PSF = padarray(PSF, [imSize-ap_field, imSize-ap_field]/2, 0, 'both');
PSF = PSF./sum(PSF(:));

stim = fftshift(ifft2(fft2(stimTemplate).*fft2(PSF)));
