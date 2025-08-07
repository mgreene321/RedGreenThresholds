% for each location, calculate euclidean distance from every location
% (including itself)

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);

for subj = 1:numel(subjects)
    [xm{subj}, ym{subj}] = findTargetedLocationsInMaster(subjects{subj});
    allLocs{subj} = [transpose(xm{subj}) transpose(ym{subj})];
    for i = 1:size(xm{subj},2)
        loc = allLocs{subj}(i,:);
        distanceMatrix{subj}(:,i) = sqrt(sum((loc - allLocs{subj}).^2, 2));



    end


end