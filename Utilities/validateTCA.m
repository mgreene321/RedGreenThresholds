% Go through all found crosses on good trials, ensure that red or green
% cross locs + tca offsets get you to IR cross (or close to it
filepath = 'Z:\Local_Share\Imaging_Data\AOVIS\';
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
RawDataTable = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);


counter = 1;

appliedTCA = [];
measuredTCA = [];

%stimOffs = [];
for i = 1:size(ADT,1)
    subjectId = ADT.SubjectID{i};
    expFolder = ADT.Folder{i};
    expDate = regexp(expFolder, '\d{1,2}_\d{1,2}_\d{4}', 'match');
    expTimes = ADT.ExperimentTimes{i};


    for T = 1:numel(expTimes)
        expTime = expTimes{T};
        % load good delivery data
        goodDeliveries = importdata(fullfile(RedGreenThresholdsPath,  subjectId, expFolder, expTime, 'goodDeliveries.mat'));
        goodTrialNumbers = unique(goodDeliveries(:,4,1));
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

        % figure out tca
        dataDir = dir([expTimeFolder, '\RedGreenThresholds*.mat']);
        load(fullfile(dataDir.folder, dataDir.name));
        data = exp_data.data_matrix;

        tca = [exp_data.expParams.redTCA_X...
            exp_data.expParams.redTCA_Y...
            exp_data.expParams.greenTCA_X...
            exp_data.expParams.greenTCA_Y];

        [xx yy] = meshgrid(0:exp_data.expParams.gridHeight - 1, 0:exp_data.expParams.gridWidth - 1);
        stimOff = max(exp_data.expParams.stimDiam) .* [xx(:) yy(:)];

        % go through each trial


%         RawDataTableFragment = RawDataTable(strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.Folder, expDateFolder{1})...
% & strcmpi(RawDataTable.ExperimentTime, expTime),:);


        for j = 1:numel(goodTrialNumbers)
            aom = data(goodTrialNumbers(j),3);
            loc = data(goodTrialNumbers(j),1);
            trialNum = goodTrialNumbers(j);
            IRcross = goodDeliveries(goodDeliveries(:,4) == trialNum, 1:2, 1);
            redCross = goodDeliveries(goodDeliveries(:,4) == trialNum, 1:2, 2);
            greenCross = goodDeliveries(goodDeliveries(:,4) == trialNum, 1:2, 3);

   
            % temp_appliedRedTCA = redCross - IRcross;
            % temp_appliedGreenTCA = greenCross - IRcross;

            if aom == 1
                temp_appliedTCA = redCross - IRcross - stimOff(loc,:);
                appliedTCA = [appliedTCA; temp_appliedTCA];
                temp_measuredTCA = repmat(tca(1:2), size(temp_appliedTCA,1), 1);
            elseif aom == 2
                temp_appliedTCA = greenCross - IRcross - stimOff(loc,:);
                appliedTCA = [appliedTCA; temp_appliedTCA];
                temp_measuredTCA = repmat(tca(3:4), size(temp_appliedTCA,1), 1);
            end
            measuredTCA = [measuredTCA; temp_measuredTCA];
            % 
            % appliedTCA = [appliedTCA; temp_appliedTCA];
            % measuredTCA = [measuredTCA; temp_measuredTCA];

          allExpDates{counter,:} = expDate{:};
             counter = counter + 1;


        end


        % for j = 1:size(goodDeliveries,1)
        % 
        % 
        %     IRcross = goodDeliveries(j, 1:2, 1);
        %     redCross = goodDeliveries(j,1:2,2);
        %     greenCross = goodDeliveries(j,1:2,3);
        % 
        %     appliedTCA(counter,1:2) = redCross - IRcross;
        %     appliedTCA(counter,3:4) = greenCross - IRcross;
        % 
        %     measuredTCA(counter,:) = tca;
        % 
        %        allExpDates{counter,:} = expDate{:};
        %     counter = counter + 1;
        % 
        % 
        % end




    end



end