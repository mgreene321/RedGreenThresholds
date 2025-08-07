function stimLocs = getStimLocCoordsOnMaster(subjectId, expDate, plotFlag, newFig)

if nargin < 4
newFig = 1;
elseif nargin < 3
pltotFlag = 0;
newFig = 0;

end

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
% Get OCT image and cone locs
masterData = importdata(fullfile(subjectFolder, 'OCTMaster', 'coneData.mat'));
masterConeLocs = str2double(masterData.coneData(:,[1 2]));

% Load daily data

dailyData = importdata(fullfile(subjectFolder, expDate, 'TransferFromOCT', 'dailyData.mat'));
dailyConeLocs = str2double(dailyData(:, [1 2]));
badIdx = any(isnan(dailyConeLocs),2) | any(dailyConeLocs < 0,2);

dailyCorrPts = dailyConeLocs(~badIdx,:);
masterCorrPts = masterConeLocs(~badIdx,:);


stimLocsDataDir = dir(fullfile(subjectFolder, expDate, '*\deliveryStatistics.mat'));

dailyStimLocs = [];
for i = 1:size(stimLocsDataDir,1)

temp = importdata(fullfile(stimLocsDataDir(i).folder, stimLocsDataDir(i).name));
dailyStimLocs = [dailyStimLocs; temp.meanDeliveryLocs];

end

[x_master, y_master] = transformImageCoordinates(dailyStimLocs(:,1), dailyStimLocs(:,2), 1, [], dailyCorrPts, masterCorrPts);

stimLocs = [x_master y_master];


if plotFlag
if newFig
    figure(1), title('OCT image'); hold on;
 imagesc(masterData.Raw_image_fundus_scaled), colormap gray; axis square off
    
    %add cone labels
    for i = 1:size(masterConeLocs,1)
        if masterData.coneData(i,3) == "L"
            plot(masterConeLocs(i,1), masterConeLocs(i,2), 'r.', 'MarkerSize', 9)
        elseif masterData.coneData(i,3) == "M"
            plot(masterConeLocs(i,1), masterConeLocs(i,2), 'g.', 'MarkerSize', 9)
        elseif masterData.coneData(i,3) == "S"
            plot(masterConeLocs(i,1), masterConeLocs(i,2), 'b.', 'MarkerSize', 9)
        else
            plot(masterConeLocs(i,1),masterConeLocs(i,2), 'k.', 'MarkerSize', 9)
        end
    end

else
end
  
figure(1);
scatter(x_master, y_master, 'yo'), set(gca, 'ydir', 'reverse')

figure;  title('Daily Image'); hold on;
dailyImage = importdata(fullfile(subjectFolder, expDate, 'TransferFromOCT', 'meanFrame.mat'));
imagesc(dailyImage), colormap gray, axis square off

set(gca, 'ydir', 'reverse')

    for i = 1:size(dailyConeLocs,1)
        if dailyData(i,3) == "L"
            plot(dailyConeLocs(i,1), dailyConeLocs(i,2), 'r.', 'MarkerSize', 9)
        elseif dailyData(i,3) == "M"
            plot(dailyConeLocs(i,1), dailyConeLocs(i,2), 'g.', 'MarkerSize', 9)
        elseif dailyData(i,3) == "S"
            plot(dailyConeLocs(i,1), dailyConeLocs(i,2), 'b.', 'MarkerSize', 9)
        else
            plot(dailyConeLocs(i,1),dailyConeLocs(i,2), 'k.', 'MarkerSize', 9)
        end
    end

     missingConeLocs = importdata(fullfile(subjectFolder, expDate, 'TransferFromOCT', 'likelyMissingCones.mat'));
     plot(missingConeLocs(:,1), missingConeLocs(:,2), 'k.', 'MarkerSize', 9);

scatter(dailyStimLocs(:,1), dailyStimLocs(:,2), 'yo');


end
end