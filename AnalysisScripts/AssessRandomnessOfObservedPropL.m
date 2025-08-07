seed = 1;
rng(seed);
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
analyzedDataTable = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);

%% Create PSF
pupilDiamMm = 6.5;
ap_field = 512;
psf_pupil = pupilDiamMm;
zernike_pupil = pupilDiamMm;
field_size = 0.9*60;
diff_limited = 1;
defocus = 0.05;
lambdaUm = 543/1000;
PSF543 = generate_PSF(ap_field, psf_pupil, zernike_pupil, field_size, lambdaUm, diff_limited, defocus);
PSF543 = PSF543./sum(PSF543(:)); % normalize volume to 1
PSF543 = padarray(PSF543, [100 100], 0, 'both');

%% Convolve stimulus with PSF
imSize = 712;
stimDiamPix = 21;
stimTemplate = zeros(imSize, imSize);
stimTemplate((imSize/2 + 1) - fix(stimDiamPix/2): (imSize/2 + 1) + fix(stimDiamPix/2), (imSize/2 + 1) - fix(stimDiamPix/2):(imSize/2 + 1) + fix(stimDiamPix/2)) = 1;
stim543 = fftshift(ifft2(fft2(stimTemplate).*fft2(PSF543)));

%% Compute cone aperture function
coneSpacingPixels = 11;
coneDiamPixels = coneSpacingPixels;
coneAperturePixels = coneDiamPixels.*.5;
gaussianSigma = coneAperturePixels./2.355;% = 0.2502 arcmin
windowSize = imSize;
coneAperture = fspecial('gaussian', windowSize, gaussianSigma);

subjects = unique(analyzedDataTable.SubjectID);

% Go through for each subject and run simulation, using latest
% coneApertures mat file

for subj = 1:numel(subjects)
    N = []; LLMsamples = [];
    %find latest coneApertures.mat
    subjectSubFolders = {dir(fullfile(RedGreenThresholdsPath, subjects{subj})).name}';
    dateFolderIdx = cellfun(@(x) ~isempty(x),regexp(subjectSubFolders, '^\d{1,2}_\d{1,2}_\d{4}'));
    datesFolders = subjectSubFolders(dateFolderIdx);
    datesStr = regexp(subjectSubFolders, '^\d{1,2}_\d{1,2}_\d{4}', 'Match');
    datesStr = [datesStr{:}]';

    dates = cellfun(@(x) datetime(x, 'InputFormat', 'MM_dd_yyyy'), datesStr);
    latestDateIdx = dates == max(dates);
    latestDateFolder = datesFolders(latestDateIdx);

    coneAperturesFile = fullfile(RedGreenThresholdsPath, subjects{subj}, latestDateFolder, 'TransferFromOCT', 'coneApertures.mat');
    coneApertures = importdata(coneAperturesFile{1});

    L = coneApertures(:,:,1);
    M = coneApertures(:,:,2);
    Lconv = fftshift(ifft2(fft2(stim543).*fft2(L)));
    Mconv = fftshift(ifft2(fft2(stim543).*fft2(M)));
    LLM = Lconv./(Lconv+Mconv);
    LLM_crop = imcrop(LLM, [100 100 412 412]);
    LLM_vec = LLM_crop(:);


    % sample as many locations as were actually targeted, 10,000 times

    numTargetedLocs = sum(strcmpi(analyzedDataTable.SubjectID, subjects{subj}));
    for i = 1:10000
        qq = randi([1 numel(LLM_crop)], 1, numTargetedLocs);
        LLMsamples(i,:) = LLM_vec(qq);
    end
    binEdges = 0:(1/8):1;
    midPoints = binEdges(1:end-1) + diff(binEdges)/2;

    for i = 1:10000
        N(i,:) = histcounts(LLMsamples(i,:), binEdges);
    end

    Nmean = mean(N,1);
    Nstd = std(N,1);
    Nsem = Nstd./sqrt(1e4);

    propL = analyzedDataTable.L./(analyzedDataTable.L + analyzedDataTable.M);
    figure; hold on, title(subjects{subj}); xlabel('Proportion L'); ylabel('Count');

    Nmean2 = repelem(Nmean,2);
    Nstd2 = repelem(Nstd,2);
    newX = [binEdges(1) repelem(binEdges(2:end-1),2) binEdges(end)];
    x = [newX, fliplr(newX)];
    y1 = Nmean2 + Nstd2;
    y2 = Nmean2 - Nstd2;
    inBetween = [y1, fliplr(y2)];
    fill(x, inBetween, 'r', 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'LineWidth', 1.5);
    plot(newX, y1, 'LineStyle', 'none', 'Color', 'r', 'LineWidth', 1.5);

    %plot(newX, Nmean2, 'Color', 'r', 'LineWidth', 3)
    for i = 1:numel(newX)-1
        if newX(i) ~= newX(i+1)
            plot([newX(i) newX(i+1)], [Nmean2(i) Nmean2(i+1)], 'Color', 'r', 'LineWidth', 3)
        else
        end
    end


    plot(newX, y2, 'LineStyle', 'none', 'Color', 'r', 'LineWidth', 1.5)

    histHand(subj) =histogram(propL(strcmpi(analyzedDataTable.SubjectID, subjects{subj})), 'BinEdges', binEdges);
    histHand(subj).LineWidth = 1.5;
    histHand(subj).FaceColor ='none';

    %errorbar(midPoints, Nmean, Nstd, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')

    ylim([0 13])
    histAxes(subj) = gca;

    % compare simulated and empirical distribution
    


end