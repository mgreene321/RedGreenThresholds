function dailyData = getDailyConeData(subjectId, expDate)
addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\AOMcontrol\'));
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);

load(fullfile(subjectFolder, 'OCTMaster', 'coneData.mat'));
load(fullfile(subjectFolder, expDate, 'TransferFromOCT', 'meanFrame.mat'));


master_cone_locs = str2double(coneData(:, [1 2]));
cone_labels = coneData(:,3);
x_master = master_cone_locs(:,1);
y_master = master_cone_locs(:,2);

[x_daily, y_daily] = transformImageCoordinates(x_master,y_master, 0,0);
daily_cone_locs = [x_daily y_daily];
initialDailyData = [daily_cone_locs coneData(:,3)];



% Manually tweak cone locations, especially those around where the stimulus
% was delivered

%Get all the stimulus delivery locations

tempDir = dir(fullfile(subjectFolder, expDate, '**\deliveryStatistics.mat'));
allMeanDeliveryLocs = [];

for i = 1:size(tempDir,1)
    load(fullfile(tempDir(i).folder, tempDir(i).name), 'meanDeliveryLocs');
    allMeanDeliveryLocs = [allMeanDeliveryLocs; meanDeliveryLocs];
end


daily_cone_locs = editConeLocations(meanFrame, daily_cone_locs, cone_labels, allMeanDeliveryLocs);
cone_centers = centerCones(meanFrame, daily_cone_locs);

dailyData = [cone_centers coneData(:,3)];
meanFrame(isnan(meanFrame)) = 0;
[allConesInClassRegion, isLabeled] = identifyUnlabeledCones(meanFrame, cone_centers);
likelyMissingCones = allConesInClassRegion(~isLabeled,:);


likelyMissingCones = centerCones(meanFrame, likelyMissingCones);

saveFolder = fullfile(subjectFolder, expDate, 'TransferFromOCT');

save(fullfile(saveFolder, 'initialDailyData.mat'), 'initialDailyData');
save(fullfile(saveFolder, 'dailyData.mat'), 'dailyData');
save(fullfile(saveFolder, 'likelyMissingCones.mat'), 'likelyMissingCones');
end