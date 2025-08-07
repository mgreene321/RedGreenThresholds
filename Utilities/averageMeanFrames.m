function meanFrame = averageMeanFrames(subjectId, expDate)

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
tempDir = dir(fullfile(subjectFolder,expDate, '**\findAllCrossesOutput.mat'));

TransferFromOCTFolder = 'TransferFromOCT';

for i = 1:size(tempDir,1)
    load(fullfile(tempDir(i).folder, 'findAllCrossesOutput.mat'))
    meanFrames(:,:,i) = meanFrame;

end

if ~isdir(fullfile(subjectFolder, expDate, TransferFromOCTFolder))
    mkdir(fullfile(subjectFolder, expDate, TransferFromOCTFolder));
else
end

meanFrame = mean(meanFrames, 3, 'omitnan');
imwrite(meanFrame./max(meanFrame(:)), fullfile(subjectFolder, expDate, TransferFromOCTFolder, 'meanFrame.tif'));
save(fullfile(subjectFolder, expDate,TransferFromOCTFolder, 'meanFrame.mat'), 'meanFrame'); 