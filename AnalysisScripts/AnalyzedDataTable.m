% Generates table with the following variables:

%   -SubjectID
%   -Folder: named for experiment day in the format mm_dd_yyyy
%   -ExperimentTimes: cell array containing experimental block times in the
%   format HH_MM_SS
%   -Block: index of experimental block, starting with 1
%   -Location: index of retinal location within block
%   -StimDiamPix: stimulus diameter in pixels
%   -StimDurFrames: stimulus duraion in frames
%   -RedThresholdPowerAU: 680 nm threshold in 0-1 arbitrary intensity units
%   -GreenThresholdPowerAU: 543 nm threshold in 0-1 arbitrary intensity
%   units
%   -RedPowerWattsAt680nm: maximum power [W] measured at 680 nm
%   -GreenPowerWattsAt543nm: maximum power [W] measured at 543 nm
%   -RedSlope: slope of fitted psychometric function for 680 nm
%   -GreenSlope: slope of fitted psychometric function for 543 nm
%   -L: mean number of stimulated L-cones for targeted location
%   -M: mean number of stimulated M-cones for targeted location
%   -S: mean number of stimulated S-cones for targeted location
%   -NumRed543: number of red responses to 543 nm stimulus
%   -NumRed680: number of red responses to 680 nm stimulus
%   -NumGreen543: number of green responses to 543 nm stimulus
%   -NumGreen680: number of green resposnes to 680 nm stimulus
%   -NumAchrom543: number of achromatic responses to 543 nm stimulus
%   -NumAchrom680: number of achromatic responses to 680 nm stimulus
%   -NumNotSeen543: number of not seen responses to 543 nm stimulus
%   -NumNotSeen680 number of not seen responses to 680 nm stimulus

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
% Load raw data

%RawDataTable = RedGreenThresholds_DataTable;
RawDataTable = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\RawDataTable.mat');

subjects = unique(RawDataTable.SubjectID);
stimDiams = unique(RawDataTable.StimDiamPix);
stimDurs = unique(RawDataTable.StimDurFrames);

% For each subject, get analysis outputs

analysisOutputDirs = [];
for subj = 1:size(subjects,1)
    subjIdx = strcmpi(RawDataTable.SubjectID, subjects{subj});
    expDateFolders{subj} = unique(RawDataTable.Folder(subjIdx));
    
    for D = 1:size(expDateFolders{subj},1)
        temp = dir(fullfile(RedGreenThresholdsPath, subjects{subj},expDateFolders{subj}{D}, 'analysisOutput*.mat'));
        analysisOutputDirs{subj, D} = temp;
    end
end

% Now for each subject and date, load the analysis output dir, find the
% expTimes for the analysis, get the corresponding data from the raw table
% and average L, M, S, and color responses.
% Subject Date expTime Location L M S numRed numGreen numAchrom numNotSeen

analyzedDataTable = [];

% for each subject
for subj = 1:size(subjects,1)
    subjIdx = strcmpi(RawDataTable.SubjectID, subjects{subj});
    
    %for each experiment date
    for D = 1:size(expDateFolders{subj},1)
        Folder = expDateFolders{subj}(D);
        % for each analysisOutput.mat file
        for i = 1:size(analysisOutputDirs{subj, D},1)
            Block = i;
            % load analysisOutput.mat
            analysisOutput = importdata(fullfile(analysisOutputDirs{subj, D}(i).folder, analysisOutputDirs{subj, D}(i).name));
            ExperimentTimes = analysisOutput.expTime;
            
            % Start setting up table
            
            SubjectID = subjects(subj);
            
            
            for diam = 1:numel(stimDiams)
                StimDiamPix = stimDiams(diam);
                for dur = 1:numel(stimDurs)
                    StimDurFrames = stimDurs(dur);
                    tempTbl = RawDataTable(strcmpi(RawDataTable.SubjectID, subjects{subj})...
                        & strcmpi(RawDataTable.Folder, Folder)...
                        & ismember(RawDataTable.ExperimentTime, ExperimentTimes)...
                        & RawDataTable.StimDiamPix == StimDiamPix...
                        & RawDataTable.StimDurFrames == StimDurFrames,:);
                    
                    
                    % For each location in tempTbl, compute color resposes
                    
                    for loc = 1:max(tempTbl.Location)
                        Location = loc;
                        % RedThreshold_au = 10.^(analysisOutput.psiParamsFit{loc,1}(1)./20);
                        % GreenThreshold_au = 10.^(analysisOutput.psiParamsFit{loc,2}(1)./20);
                        
                        % Delivery stats
                        
                        meanXc = mean(tempTbl.Xc(tempTbl.Location == loc));
                        meanYc = mean(tempTbl.Yc(tempTbl.Location == loc));
                        
                        stdXc = std(tempTbl.Xc(tempTbl.Location == loc));
                        stdYc = std(tempTbl.Yc(tempTbl.Location == loc));
                        
                        try
                        RedThresholdPowerAU = analysisOutput.redThresholdPowerAU{diam,dur}(loc);
                        GreenThresholdPowerAU = analysisOutput.greenThresholdPowerAU{diam,dur}(loc);
                        catch
                             RedThresholdPowerAU = analysisOutput.redThresholdPowerAU(loc);
                        GreenThresholdPowerAU = analysisOutput.greenThresholdPowerAU(loc);
                        end
                        
                        RedPowerWattsAt680nm = analysisOutput.redPowerWattsAt680nm;
                        GreenPowerWattsAt543nm = analysisOutput.greenPowerWattsAt543nm;
                        
                        try
                        RedSlope = analysisOutput.psiParamsFit{loc,1, diam, dur}(2);
                        GreenSlope = analysisOutput.psiParamsFit{loc,2, diam, dur}(2);
                        catch
                            RedSlope = analysisOutput.psiParamsFit{loc,1}(2);
                        GreenSlope = analysisOutput.psiParamsFit{loc,2}(2);
                        end
                        
                        %                 RedThreshold_J = analysisOutput.thresholds(loc,1);
                        %                 GreenThreshold_J = analysisOutput.thresholds(loc,2);
                        
                        NumRed680 = sum(tempTbl.Color(tempTbl.Location == loc) == 1 & tempTbl.Channel(tempTbl.Location == loc) == 1);
                        NumRed543 = sum(tempTbl.Color(tempTbl.Location == loc) == 1 & tempTbl.Channel(tempTbl.Location == loc) == 2);
                        NumGreen680 = sum(tempTbl.Color(tempTbl.Location == loc) == 2 & tempTbl.Channel(tempTbl.Location == loc) == 1);
                        NumGreen543 = sum(tempTbl.Color(tempTbl.Location == loc) == 2 & tempTbl.Channel(tempTbl.Location == loc) == 2);
                        NumAchrom680 = sum(tempTbl.Color(tempTbl.Location == loc) == 3 & tempTbl.Channel(tempTbl.Location == loc) == 1);
                        NumAchrom543 = sum(tempTbl.Color(tempTbl.Location == loc) == 3 & tempTbl.Channel(tempTbl.Location == loc) == 2);
                        NumNotSeen680 = sum(tempTbl.Color(tempTbl.Location == loc) == 4 & tempTbl.Channel(tempTbl.Location == loc) == 1);
                        NumNotSeen543 = sum(tempTbl.Color(tempTbl.Location == loc) == 4 & tempTbl.Channel(tempTbl.Location == loc) == 2);
                        
                        L = mean(tempTbl.L(tempTbl.Location == loc), 'omitnan');
                        M = mean(tempTbl.M(tempTbl.Location == loc), 'omitnan');
                        S = mean(tempTbl.S(tempTbl.Location == loc), 'omitnan');
                        X = mean(tempTbl.X(tempTbl.Location == loc), 'omitnan');
                        
                        tableRow = table(SubjectID, Folder, {ExperimentTimes}, Block, Location,...
                            StimDiamPix, StimDurFrames,...
                            RedThresholdPowerAU, GreenThresholdPowerAU, RedPowerWattsAt680nm, GreenPowerWattsAt543nm,...
                            meanXc, meanYc, stdXc, stdYc, RedSlope, GreenSlope,...
                            L, M, S,X,...
                            NumRed543, NumRed680, NumGreen543, NumGreen680, NumAchrom543, NumAchrom680, NumNotSeen543, NumNotSeen680);
                        tableRow.Properties.VariableNames{3} = 'ExperimentTimes';
                        analyzedDataTable = [analyzedDataTable; tableRow];
                    end
                end
                
            end
        end
    end
end