function indices = getLocationIndicesInRawDataTable
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
%Load analyzed data table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
% Load in raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);
indices = cell(size(ADT,1),1);
for i = 1:size(ADT,1)
    subjectId = ADT.SubjectID{i};
    expFolder = ADT.Folder{i};
    expTimes = ADT.ExperimentTimes{i};
    location = ADT.Location(i);


    indices{i} = strcmpi(subjectId, RDT.SubjectID) &...
                 strcmpi(expFolder, RDT.Folder) &...
                 ismember(RDT.ExperimentTime, expTimes) &...
                 RDT.Location == location;
end