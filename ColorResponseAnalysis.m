%analyzeNearThresholdData
vidPath = 'Z:\Local_Share\Imaging_Data\AOVIS\';

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);

RawDataTable = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
AnalyzedDataTable;
subjects = unique(RawDataTable.SubjectID);

for subj = 1:numel(subjects)
%figure(subj); hold on
    subjectId = subjects{subj};
    expFolders = unique(RawDataTable.Folder(strcmpi(RawDataTable.SubjectID, subjects{subj})));
    
    locCounter = 1;
    for F = 1:numel(expFolders)
        
        expFolder = expFolders{F};
        
        Blocks = unique(analyzedDataTable.Block(strcmpi(analyzedDataTable.SubjectID, subjectId) & strcmpi(analyzedDataTable.Folder, expFolder),:));
        
        
        for B = 1:numel(Blocks)
            % get data from raw data table and from analyzed data
            
            % for each location
            % divide stimulus intensities by threshold estimate
            % plot cumulative "Red", "Green", "Achromatic", "Not Seen" responses
            % vs. intensity in threshold units
            
            Block = Blocks(B);
            idx1 = strcmpi(analyzedDataTable.SubjectID, subjectId) & strcmpi(analyzedDataTable.Folder, expFolder) & analyzedDataTable.Block == Block;
            analyzedData = analyzedDataTable(idx1,:);
            
            expTimes = unique([analyzedData.ExperimentTimes{:}]);
            
            redThresholds = analyzedData.RedThreshold_au;
            greenThresholds = analyzedData.GreenThreshold_au;
            
%             redAuToJ = analyzedData.RedThreshold_J./analyzedData.RedThreshold_au;
%             greenAuToJ = analyzedData.GreenThreshold_J./analyzedData.GreenThreshold_au;
            
            L = analyzedData.L; M = analyzedData.M; S = analyzedData.S;
            propL = L./(L+M);
            
            minIntensity = 10^-0.5;
            maxIntensity = 10^0.5;
            for loc = 1:max(analyzedData.Location)
                % get trial by trial data
                idx2 = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.Folder, expFolder) & ismember(RawDataTable.ExperimentTime, expTimes) & RawDataTable.Location == loc;
                tempRawData = RawDataTable(idx2,:);
                
                % convert intensities to threshold units
                tempRawData.Intensity(tempRawData.Channel == 1) = tempRawData.Intensity(tempRawData.Channel == 1)./ redThresholds(loc);
                tempRawData.Intensity(tempRawData.Channel == 2) = tempRawData.Intensity(tempRawData.Channel == 2) ./  greenThresholds(loc);
                
                % Compute P("red" | 680) for near threshold stimuli
                
                Red680 = tempRawData.Channel == 1 & tempRawData.Color == 1 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                Green680 = tempRawData.Channel == 1 & tempRawData.Color == 2 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                Achrom680 = tempRawData.Channel == 1 & tempRawData.Color == 3 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                NS680 = tempRawData.Channel == 1 & tempRawData.Color == 4 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                
                P_Red680{subj}(locCounter,:) = sum(Red680) ./ (sum(Green680) + sum(Achrom680) + sum(NS680));
                
                %P_Red680 = sum(Red680);
                % Compute P("green" | 543)
                
                Green543 = tempRawData.Channel == 2 & tempRawData.Color == 2 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                Red543 = tempRawData.Channel == 2 & tempRawData.Color == 1 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                Achrom543 = tempRawData.Channel == 2 & tempRawData.Color == 3 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                NS543 = tempRawData.Channel == 2 & tempRawData.Color == 4 & (tempRawData.Intensity >= minIntensity & tempRawData.Intensity <= maxIntensity);
                
                P_Green543{subj}(locCounter,:) = sum(Green543) ./ (sum(Red543) + sum(Achrom543) + sum(NS543));
                %P_Green680 = sum(Green543);
                
                locCounter = locCounter + 1;
                
            end
          %  scatter(P_Green543, P_Red680, 'filled');
%             P_Red680 = [];
%             P_Green543 = [];
        end
    end
end
