function numTrials = getTotalTrialsPerLocation

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
analyzedDataTable = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = analyzedDataTable(analyzedDataTable.StimDurFrames == 3,:);
numTrials = nan(size(ADT,1),1);

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

        stim_matrix =tempData.exp_data.stim_matrix;
        if size(stim_matrix,2) == 6
            stimDiamIdx = find(tempData.exp_data.expParams.stimDiam == 21);
            stimDurIdx = find(tempData.exp_data.expParams.stimDur == 3);
            stim_matrix = stim_matrix(stim_matrix(:,4) == stimDiamIdx & stim_matrix(:,5) == stimDurIdx,:);


        end


        N = sum(stim_matrix(:,1) == ADT.Location(i));
        tempN = tempN + N;
    end
numTrials(i,:) = tempN;    
end


end