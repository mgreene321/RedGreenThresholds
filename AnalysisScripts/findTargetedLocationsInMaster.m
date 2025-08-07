function [xm_mean, ym_mean] = findTargetedLocationsInMaster(subjectId)

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
RawDataTable = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'AnalyzedDataTable.mat'));

analyzedDataTable = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);
subjectRawData= RawDataTable(strcmpi(RawDataTable.SubjectID, subjectId),:);
subjectAnalyzedData = analyzedDataTable(strcmpi(analyzedDataTable.SubjectID, subjectId),:);

OCTdata = importdata(fullfile(RedGreenThresholdsPath, subjectId, 'OCTMaster', 'coneData.mat'));
coneData = OCTdata.coneData;

OCTimage = OCTdata.Raw_image_fundus_scaled;
figure, imshow(OCTimage),hold on
%set(gca, 'Units', 'pixels');
%set(gca, 'Position', [1 1 size(OCTimage,1) size(OCTimage,2)]);
%set(gcf, 'Units', 'pixels');
%set(gcf, 'Position', [1 1 round(1.1*size(OCTimage,1)) round(1.1*size(OCTimage,2))]);



% Plot cone markers
master_cone_locs = str2double(coneData(:,1:2));
cone_labels = coneData(:,3);
plot(master_cone_locs(cone_labels == "L",1), master_cone_locs(cone_labels == "L",2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerSize', 5);
plot(master_cone_locs(cone_labels == "M",1), master_cone_locs(cone_labels == "M",2),'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.7 0.1], 'MarkerSize', 5);
plot(master_cone_locs(cone_labels == "S",1), master_cone_locs(cone_labels == "S",2), 'LineStyle', 'none','Marker', 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.2 0.2 1], 'MarkerSize', 5);


% For each date folder
dateFolders = unique(subjectAnalyzedData.Folder);

ii = 1;
for F = 1:size(dateFolders,1)
    dateFolderIdx = strcmpi(subjectAnalyzedData.Folder, dateFolders{F});
    
    % load cone transfer data
    dailyData = importdata(fullfile(RedGreenThresholdsPath, subjectId, dateFolders{F}, 'TransferFromOCT', 'dailyData.mat'));
    
    dailyConeLocs = str2double(dailyData(:, [1 2]));
    badIdx = any(isnan(dailyConeLocs),2) | any(dailyConeLocs < 0,2);
    dailyCorrPts = dailyConeLocs(~badIdx,:);
    masterConeLocs = str2double(OCTdata.coneData(:, [1 2]));
    masterCorrPts = masterConeLocs(~badIdx,:);
    % For each block
    for B = 1:max(subjectAnalyzedData.Block(dateFolderIdx))
        blockIdx = dateFolderIdx & subjectAnalyzedData.Block == B;
        temp = subjectAnalyzedData(blockIdx,:);
        for loc = 1:max(subjectAnalyzedData.Location(blockIdx))
            % % Get sensitivity ratio
            % redSensitivity = 1./(subjectAnalyzedData.RedThreshold_J(blockIdx & subjectAnalyzedData.Location == loc));
            % greenSensitivity = 1./(subjectAnalyzedData.GreenThreshold_J(blockIdx & subjectAnalyzedData.Location == loc));
            %
            %  RGratio= redSensitivity./greenSensitivity;
            
            idx = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.Folder, dateFolders{F})...
                & ismember(RawDataTable.ExperimentTime, unique([temp.ExperimentTimes{:}])) & RawDataTable.Location == loc;
            
            stimX = RawDataTable.Xc(idx); stimY = RawDataTable.Yc(idx);
            [xm, ym] = transformImageCoordinates(stimX, stimY, 1, [], dailyCorrPts, masterCorrPts);
            
            xm_mean(ii) = mean(xm); ym_mean(ii) = mean(ym);
            
            %rectColor = [(RGratio - minRGratio)./(maxRGratio - minRGratio) 1-(RGratio - minRGratio)./(maxRGratio - minRGratio) 0];
            rectColor = 'w';
            drawrectangle('Position', [xm_mean(ii)-11 ym_mean(ii)-11 21 21], 'LineWidth', 1.5, 'Color', rectColor, 'FaceAlpha', 0.3, 'InteractionsAllowed', 'none');
            
            ii = ii + 1;
        end
    end
    
end

xlim([0 600]); ylim([0 600]); 
axis equal
end