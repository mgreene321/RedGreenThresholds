function RawDataTable = RedGreenThresholds_DataTable
vidPath = 'Z:\Local_Share\Imaging_Data\AOVIS\';

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);

% SubjectId expDate expTime Location StairNum AOM StimIntensity YesNo Color GoodTrial


RedGreenThresholdsDir = dir(RedGreenThresholdsPath);
folderNames = {RedGreenThresholdsDir(:).name}';
temp = regexp(folderNames, '\d{5}[LR]');
subjects = folderNames(cellfun(@(x) ~isempty(x), temp));

% GET RID OF HANNAH
subjects(strcmpi(subjects, '20225R')) = [];
subjects(strcmpi(subjects, '10003L')) = [];

%for each subject
RG_tbl = [];
for subj = 1:size(subjects,1)
    subjectId = subjects{subj};
    subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
    
    % Find dates
    expDateDir = dir(subjectFolder);
    folderNames = {expDateDir(:).name}';
    temp = regexp(folderNames, '^\d{1,2}_\d{1,2}_\d{4}');
    expDateFolders =  folderNames(cellfun(@(x) ~isempty(x), temp));
    expDates = regexp(expDateFolders, '^\d{1,2}_\d{1,2}_\d{4}', 'Match');
    expDates = [expDates{:}]';
    
    % For eahc date, find times
    for D = 1:size(expDateFolders,1)
        expDate = expDates{D};
        expDateFolder = fullfile(subjectFolder,expDateFolders{D});
        expTimeDir = dir(expDateFolder);
        folderNames = {expTimeDir(:).name}';
        temp = regexp(folderNames, '\d{1,2}_\d{1,2}_\d{1,2}');
        expTimes = folderNames(cellfun(@(x) ~isempty(x), temp));
        
        for T = 1:size(expTimes,1)
            expTime = expTimes{T};
            expTimeFolder = fullfile(expDateFolder, expTime);
            % load in experiment data
            expVideoFolder = fullfile(vidPath, subjectId, [expDate '_' expTime]);
            
            if ~isdir(expVideoFolder)
                temp = regexp(subjectId, '\d{5}|[LR]', 'match');
                subjectNumber = temp{1};
                subjectEye = temp{2};
                if strcmpi(subjectEye, 'L')
                    expVideoFolder = fullfile(vidPath, subjectNumber, 'Left', [expDate '_' expTime]);
                else
                    expVideoFolder = fullfile(vidPath, subjectNumber, 'Right', [expDate '_' expTime]);
                end
            else
            end
            
            % Add response data to table
            expVideoDir = dir(fullfile(expVideoFolder, 'RedGreenThresholds*.mat'));
            expData = importdata(fullfile(expVideoDir(1).folder, expVideoDir(1).name));
            expData = expData.exp_data;
            
            
            if size(expData.data_matrix,2) == 6
                tempMat = expData.data_matrix(:, 1:3);
                tempMat(:, end+1) = 1;
                tempMat(:, end+1) = 1;
                tempMat(:, end+1:end+3) = expData.data_matrix(:, end-2:end);
                expData.data_matrix = tempMat;
                
            else
            end
            
            % go through and replace stim size and dur indices with the
            % actual size and durations
            
            numDiams = numel(expData.expParams.stimDiam);
            numDurs = numel(expData.expParams.stimDur);
            
            for i = 1:numDiams
                expData.data_matrix(expData.data_matrix(:,4) == i,4) = expData.expParams.stimDiam(i);
            end
            
             for i = 1:numDurs
                expData.data_matrix(expData.data_matrix(:,5) == i,5) = expData.expParams.stimDur(i);
            end

            tca = [expData.expParams.redTCA_X...
                expData.expParams.redTCA_Y...
                expData.expParams.greenTCA_X...
                expData.expParams.greenTCA_Y];

            % figure out stim offsets
       

            tempTbl1 = array2table(expData.data_matrix);
            
            tempTbl1.Properties.VariableNames = {'Location', 'StairNum', 'Channel', 'StimDiamPix', 'StimDurFrames', 'IntensityAU', 'YesNo', 'Color'};
            
            subjectIdCol = repmat({subjectId}, [size(tempTbl1,1) 1]);
            expDateCol = repmat({expDate}, [size(tempTbl1,1),1]);
            expDateFolderCol = repmat({expDateFolders{D}}, [size(tempTbl1, 1),1]);
            expTimeCol = repmat({expTime}, [size(tempTbl1,1),1]);
            trialNumbers = transpose(1:size(tempTbl1,1));
            
            tempTbl2 = table(subjectIdCol, expDateCol, expDateFolderCol, expTimeCol, trialNumbers, 'VariableNames', {'SubjectID', 'ExperimentDate', 'Folder', 'ExperimentTime', 'Trial'});
            tempTbl = [tempTbl2 tempTbl1];
            
            % Now add data on whether delivery was good,
            goodDeliveries = importdata(fullfile(expTimeFolder, 'goodDeliveries.mat'));
            goodTrialNumbers = unique(goodDeliveries(:,4,1));
            tempTbl = tempTbl(goodTrialNumbers,:); % HERE BAD TRIALS ARE FILTERED OUT
            % Get stimulus location
            findAllCrossesOutput = importdata(fullfile(expTimeFolder, 'findAllCrossesOutput.mat'));
            cross_XY = findAllCrossesOutput.cross_XY;
            cross_XY = cross_XY(ismember(cross_XY(:,4,1), goodTrialNumbers),:,:); % extract cross locations for good deliveries
            
            % add tco to cross locations
            cross_XY(:,1:2,2) = cross_XY(:,1:2,2) -  tca(1:2);
            cross_XY(:,1:2,3) = cross_XY(:,1:2,3) - tca(3:4);
            
            goodDeliveries(:,1:2,2) = goodDeliveries(:,1:2,2) -  tca(1:2);
            goodDeliveries(:,1:2,3) = goodDeliveries(:,1:2,3) - tca(3:4);
            % average cross locations per trial
            for t = 1:numel(goodTrialNumbers)
                loc = expData.data_matrix(goodTrialNumbers(t), 1);
                aom = expData.data_matrix(goodTrialNumbers(t),3);
                
                stimLoc(t,:) = mean(goodDeliveries(goodDeliveries(:,4,1) == goodTrialNumbers(t),1:2,aom + 1) - stimOff(loc,:), 1, 'omitnan');
            end
            Xc = stimLoc(:,1);
            Yc = stimLoc(:,2);
            
            %Cone data
            
%             try
                LMS = importdata(fullfile(expTimeFolder, 'LMS.mat'));
%                 L_per_trial = mean(LMS.L_per_frame,2, 'omitnan');
%                 M_per_trial = mean(LMS.M_per_frame,2, 'omitnan');
%                 S_per_trial = mean(LMS.S_per_frame,2, 'omitnan');

                  L_per_trial = LMS.L_per_trial;
                  M_per_trial = LMS.M_per_trial;
                  S_per_trial = LMS.S_per_trial;
                  X_per_trial = LMS.X_per_trial;
                L = []; M = []; S = []; X = []; stimLoc = [];
                for i = 1:size(tempTbl,1) % for each good trial
                    
%                     L(i,:) = L_per_trial(tempTbl.Location(i),:,i);
%                     M(i,:) = M_per_trial(tempTbl.Location(i),:,i);
%                     S(i,:) = S_per_trial(tempTbl.Location(i),:,i);
                      L(i,:) = L_per_trial{tempTbl.Location(i)}(i);
                      M(i,:) = M_per_trial{tempTbl.Location(i)}(i);
                      S(i,:) = S_per_trial{tempTbl.Location(i)}(i);
                      X(i,:) = X_per_trial{tempTbl.Location(i)}(i);
                end
%             catch
%                  L = []; M = []; S = []; stimLoc = [];
%                 for i = 1:size(tempTbl,1) % for each good trial
%                     
%                     L(i,:) = nan;
%                     M(i,:) = nan;
%                     S(i,:) = nan;
%                 end
%                 
%             end
%             
            
            tempTbl = [tempTbl table(L, M, S, X, Xc, Yc)];
            RG_tbl = [RG_tbl; tempTbl];
        end
        
    end
    
end
RawDataTable = RG_tbl;
save(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'), 'RawDataTable');

end