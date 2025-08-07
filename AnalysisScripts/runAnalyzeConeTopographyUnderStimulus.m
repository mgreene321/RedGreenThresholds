function runAnalyzeConeTopographyUnderStimulus

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name); % 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds'
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));

uniqueSubjects = unique(ADT.SubjectID);

for subj = 1:numel(uniqueSubjects)
    subjTable = ADT(strcmpi(ADT.SubjectID, uniqueSubjects{subj}),:);
    uniqueFolders = unique(subjTable.Folder);
    
    for fol = 1:numel(uniqueFolders)
        foldTable = subjTable(strcmpi(subjTable.Folder, uniqueFolders{fol}),:);
        uniqueExpTimes = unique(cat(1,foldTable.ExperimentTimes{:}));
        
        for et = 1:numel(uniqueExpTimes)
            expTime = uniqueExpTimes{et};
            
            analyzeConeTopographyUnderStimulus(uniqueSubjects{subj}, uniqueFolders{fol}, expTime);
            
        end
        
    end
    
    
end


end