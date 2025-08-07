function computeDeliveryError

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
% Load in raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);
 indices = getLocationIndicesInRawDataTable;

  % get ppd for each loc


 filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
for i = 1:size(ADT,1)
    expFolder = ADT.Folder{i};
    expDate = regexp(expFolder, '\d{1,2}_\d{1,2}_\d{4}', 'match');
    expTimes = ADT.ExperimentTimes{i};
    tempN = 0;
    for T = 1:numel(expTimes)
        expTime = expTimes{T};
        expTimeFolder = fullfile(filepath, ADT.SubjectID{i}, [expDate{:} '_' expTime]);

        if ~isdir(expTimeFolder)
            temp = regexp(ADT.SubjectID{i}, '\d{5}|[LR]', 'match');
            subjectNumber = temp{1};
            subjectEye = temp{2};

            if strcmpi(subjectEye, 'L')
                expTimeFolder = fullfile(filepath, subjectNumber, 'Left', [expDate{:} '_' expTime]);
            else
                expTimeFolder = fullfile(filepath, subjectNumber, 'Right', [expDate{:} '_' expTime]);
            end

        else
        end


        dataDir = dir(fullfile(expTimeFolder, 'RedGreenThresholds*.mat'));
        tempData = importdata(fullfile(dataDir.folder, dataDir.name)); 

        PPD_X(i,:) = tempData.exp_data.sessionParams.PPD_X;
        PPD_Y(i,:) = tempData.exp_data.sessionParams.PPD_Y;
 
    end
  
end

 for i = 1:size(ADT,1)

     meanXc = ADT.meanXc(i); meanYc = ADT.meanYc(i);

     xPixToArcmin = 60*(1/PPD_X(i));
     yPixToArcmin = 60*(1/PPD_Y(i));

     XcMinusMean{i} = xPixToArcmin.*(RDT.Xc(indices{i}) - meanXc); 
     YcMinusMean{i} = yPixToArcmin.*(RDT.Yc(indices{i}) - meanYc);

 end




 XcMinusMean = vertcat(XcMinusMean{:});
 YcMinusMean = vertcat(YcMinusMean{:});

 XcMinusMeanSq = XcMinusMean.^2;
 YcMinusMeanSq = YcMinusMean.^2;

 stdX = sqrt(sum(XcMinusMeanSq)./(numel(XcMinusMeanSq)-1));
 stdY = sqrt(sum(YcMinusMeanSq)./(numel(YcMinusMeanSq)-1));






end