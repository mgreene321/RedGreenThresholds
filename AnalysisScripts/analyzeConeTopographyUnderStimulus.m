function LMS = analyzeConeTopographyUnderStimulus(subjectId, expFolder, expTime)


expDate = regexp(expFolder, '\d{1,2}_\d{1,2}_\d{4}', 'Match');
expDate = expDate{:};
addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\AnalysisScripts'));
filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';

subjectNumber = char(regexp(subjectId, '\d{5}', 'match'));
subjectEye = char(regexp(subjectId, '[LR]', 'match'));

if strcmpi(subjectEye, 'L')
    subjectEye = 'Left';
else
    subjectEye = 'Right';
end

dataPath1 = fullfile(filepath, subjectId, [expDate '_' expTime]);
dataPath2 = fullfile(filepath, subjectNumber, subjectEye, [expDate '_' expTime]);

if isdir(dataPath1)
    dataPath = dataPath1;
else
    dataPath = dataPath2;
end

root = getenv('USERPROFILE');
tempDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(tempDir(1).folder, tempDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
saveFolder = fullfile(subjectFolder, expFolder, expTime);

coneTransferFolder = fullfile(subjectFolder, expFolder, 'TransferFromOCT');

% Load in experiment data

dataDir = dir([dataPath, '\RedGreenThresholds*.mat']);
load(fullfile(dataDir.folder, dataDir.name));
data = exp_data.data_matrix;

tca = [exp_data.expParams.redTCA_X...
    exp_data.expParams.redTCA_Y...
    exp_data.expParams.greenTCA_X...
    exp_data.expParams.greenTCA_Y];

% Load in cone apertures
load(fullfile(coneTransferFolder, 'coneApertures.mat'));

%coneApertures = coneApertures./max(coneApertures(:));
% Load in cross locations
load(fullfile(saveFolder, 'findAllCrossesOutput.mat'));

% Load in goodTrialNums
load(fullfile(saveFolder, 'goodDeliveries.mat'));

% For each good trial

retina_map = zeros(size(coneApertures,1,2));

% Create stimuli
% redStim = GenConvStim(680, exp_data.expParams.stimDiam, retina_map);
% greenStim = GenConvStim(543, exp_data.expParams.stimDiam, retina_map);

% PSF parameters

pupilDiamMm = 6.5;
ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05;
imSize = 712;

stimDiamPix = 21;
stimTemplate = zeros(imSize, imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;

redPSF = generate_PSF(ap_field, psf_pupil, zernike_pupil, field_size, 680/1000, diff_limited, defocus);
greenPSF = generate_PSF(ap_field, psf_pupil, zernike_pupil, field_size, 543/1000, diff_limited, defocus);

redPSF = redPSF./sum(redPSF(:));
greenPSF = greenPSF./sum(greenPSF(:));

redPSF = padarray(redPSF, [100 100], 0, 'both');
greenPSF = padarray(greenPSF, [100 100], 0, 'both');

redStim = fftshift(ifft2(fft2(stimTemplate).*fft2(redPSF)));
greenStim = fftshift(ifft2(fft2(stimTemplate).*fft2(greenPSF)));

% for each "good frame"
goodTrialNums = unique(goodDeliveries(:,4,1));

numLocs = exp_data.expParams.gridWidth .* exp_data.expParams.gridHeight;
numColors = 2;

goodData = data(goodTrialNums);
for loc = 1:numLocs     
       deliveriesPerLoc(loc, :) = sum(goodData(:,1) == loc);
end

% L = nan(numLocs, max(exp_data.expParams.stimDur), numel(goodTrialNums));
% M = L;
% S = L;

L = cell(numLocs,1);
[L{:}] = deal(nan(numel(goodTrialNums), max(exp_data.expParams.stimDur)));
M = L;
S = L;
X = L;

for trialIdx = 1:numel(goodTrialNums)
    
    % if this is a good trial
    trialNum = goodTrialNums(trialIdx);
        
          % Determine stimulus color (aom) and location index
    aom = data(trialNum,3);
    loc = data(trialNum,1);
    % Generate approximation of stimulus image on retina
        
    cross_locations = goodDeliveries(goodDeliveries(:,4,aom+1) == trialNum, 1:2, aom+1);
    
     for frm = 1:size(cross_locations,1)
        
        if aom == 1
            stim = redStim;
            tco = tca(1:2);
        elseif aom == 2
            stim = greenStim;
            tco = tca(3:4);
        end
        
      %  stim = stim./max(stim(:)); % added 6/17/24
       
%         stim = (21^2).*(stim./sum(stim(:))); % added 8/8/24
        stimLoc = cross_locations(frm,:) - tco; % changed 4/30/25
        
        % Create a BW mask that contains the PSF convolved stimulus at the
        % appropriate location
        
        xc = round(size(retina_map,2)/2); % should just be 712/2
        yc = round(size(retina_map,1)/2);
        
        %stim = circshift(stim, stimLoc - [xc yc]);
        
        try
        stim = circshift(stim, [stimLoc(2) stimLoc(1)] - [xc yc]);
        
        % Multiply mask with cone apertures
        
        L_mosaic = stim.*coneApertures(:,:,1);
        M_mosaic = stim.*coneApertures(:,:,2);
        S_mosaic = stim.*coneApertures(:,:,3);
        X_mosaic = stim.*coneApertures(:,:,4);
        
%         L(loc, frm, trialIdx) = sum(L_mosaic(:));
%         M(loc,frm, trialIdx) = sum(M_mosaic(:));
%         S(loc,frm, trialIdx) = sum(S_mosaic(:));


           L{loc}(trialIdx, frm) = sum(L_mosaic(:));
           M{loc}(trialIdx, frm) = sum(M_mosaic(:));
           S{loc}(trialIdx, frm) = sum(S_mosaic(:));
           X{loc}(trialIdx, frm) = sum(X_mosaic(:));
           
        catch
            warning('Uh oh');
        end
        
        % Take sum of 1st, 2nd, and 3rd layers of resulting matrix to get # L, M, S
        % cones
        
    end
    
end


LMS.L_per_frame = L;
LMS.M_per_frame = M;
LMS.S_per_frame = S;
LMS.X_per_frame = X;


LMS.L_per_trial = cellfun(@(x) mean(x,2, 'omitnan'), L, 'UniformOutput', false);
LMS.M_per_trial = cellfun(@(x) mean(x,2, 'omitnan'), M, 'UniformOutput', false);
LMS.S_per_trial = cellfun(@(x) mean(x,2, 'omitnan'), S, 'UniformOutput', false);
LMS.X_per_trial = cellfun(@(x) mean(x,2, 'omitnan'), X, 'UniformOutput', false);


Lnan = cellfun(@(x) all(isnan(x),2), L, 'UniformOutput', false);
Mnan = cellfun(@(x) all(isnan(x),2), M, 'UniformOutput', false);
Snan = cellfun(@(x) all(isnan(x),2), S, 'UniformOutput', false);
Xnan = cellfun(@(x) all(isnan(x),2), X, 'UniformOutput', false);

% get rid of trials where location wasn't hit (e.g., every column of row is
% nan)
for i = 1:size(L,1)
    L{i} = L{i}(~Lnan{i},:);
    M{i} = M{i}(~Mnan{i},:);
    S{i} = S{i}(~Snan{i},:);
    X{i} = X{i}(~Xnan{i},:);
    
end



% LMS.L_per_frame = L;
% LMS.M_per_frame = M;
% LMS.S_per_frame = S;


LMS.L = cellfun(@(x) mean(x,'all', 'omitnan'), L);
LMS.M = cellfun(@(x) mean(x,'all', 'omitnan'), M);
LMS.S = cellfun(@(x) mean(x,'all', 'omitnan'), S);
LMS.X = cellfun(@(x) mean(x, 'all', 'omitnan'), X);

% LMS.L = mean(sum(L,3, 'omitnan'),2, 'omitnan')./deliveriesPerLoc;
% LMS.M = mean(sum(M,3, 'omitnan'),2, 'omitnan')./deliveriesPerLoc;
% LMS.S = mean(sum(S,3, 'omitnan'),2, 'omitnan')./deliveriesPerLoc;

save(fullfile(saveFolder, 'LMS.mat'), 'LMS'); % updated 4/30/25

end