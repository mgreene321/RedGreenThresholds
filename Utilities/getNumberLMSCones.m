function [nL, nM, nS] = getNumberLMSCones(subjectId)

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);

coneDir = dir(fullfile(RedGreenThresholdsPath, subjectId, 'OCTMaster', 'Cone_identification*.mat'));
load(fullfile(coneDir.folder, coneDir.name), 'Cone_Mat_all_marked');

nL = sum(Cone_Mat_all_marked(:,3)=="L");
nM = sum(Cone_Mat_all_marked(:,3) == "M");
nS = sum(Cone_Mat_all_marked(:,3) == "S");

end