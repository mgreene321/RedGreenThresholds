root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));


ADT = ADT(ADT.StimDurFrames == 3,:);
RDT = RDT(RDT.StimDurFrames == 3,:);

folders = unique(ADT.Folder, 'stable');

coneProps10001R = [0.5545    0.3769    0.0686];
coneProps20217R = [0.7380    0.2014    0.0607 ];


% for each folder, get corresponding part of ADT and load in cone data

ii = 1;
for f = 1:size(folders,1)
    tempADT = ADT(strcmpi(folders{f}, ADT.Folder),:);
    folderName = folders{f};
    
    coneDataFile = fullfile(RedGreenThresholdsPath, unique(tempADT.SubjectID), folderName, 'TransferFromOCT', 'dailyData.mat');
    coneData = importdata(coneDataFile{1});
    
    % for each location, for each unique location run the ideal observer analysis
    uniqueBlocks = unique(tempADT.Block);
    
    for b = 1:numel(uniqueBlocks)
        block = uniqueBlocks(b);
        
        blockADT = tempADT(tempADT.Block == block,:);
        
        for loc = 1:size(blockADT,1)
            location = blockADT.Location(loc);
            
            subjectId = blockADT.SubjectID{loc};
            expTimes = blockADT.ExperimentTimes{loc};
            idx = strcmpi(RDT.SubjectID, subjectId) & RDT.Location == location & ismember(RDT.ExperimentTime, expTimes);
            
            meanXc = round(mean(RDT.Xc(idx)));
            meanYc = round(mean(RDT.Yc(idx)));
            
            stimLoc(ii,:) = [meanXc meanYc];
           
            
            if strcmpi(subjectId, '10001R')
                coneProps = coneProps10001R;
            elseif strcmpi(subjectId, '20217R')
                coneProps = coneProps20217R;
            end
            
            propL(ii,:) = blockADT.L(loc)./(blockADT.L(loc) + blockADT.M(loc));
            [thresholds543_555(ii), thresholds680_555(ii)] = IdealObserverAnalysis4(coneData, stimLoc(ii,:), coneProps);
            
            subjectIds{ii} = subjectId;
            ii = ii + 1;
            
        end
    end
    
end