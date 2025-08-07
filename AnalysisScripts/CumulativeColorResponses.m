% For given subject, experiment day, and block

RawDataTable = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\RawDataTable.mat');
AnalyzedDataTable;
subjectId = '20217R';
expFolder = '2_29_2024';
Block = 3;

% get data from raw data table and from analyzed data

% for each location
% divide stimulus intensities by threshold estimate
% plot cumulative "Red", "Green", "Achromatic", "Not Seen" responses
% vs. intensity in threshold units



idx1 = strcmpi(analyzedDataTable.SubjectID, subjectId) & strcmpi(analyzedDataTable.Folder, expFolder) & analyzedDataTable.Block == Block;
analyzedData = analyzedDataTable(idx1,:);

expTimes = unique([analyzedData.ExperimentTimes{:}]);

%redThresholds = analyzedData.RedThreshold_au;
%greenThresholds = analyzedData.GreenThreshold_au;

redAuToJ = analyzedData.RedThreshold_J./analyzedData.RedThreshold_au;
greenAuToJ = analyzedData.GreenThreshold_J./analyzedData.GreenThreshold_au;

L = analyzedData.L; M = analyzedData.M; S = analyzedData.S;
propL = L./(L+M);

for loc = 1:max(analyzedData.Location)
    % get trial by trial data
    idx2 = strcmpi(RawDataTable.SubjectID, subjectId) & strcmpi(RawDataTable.Folder, expFolder) & ismember(RawDataTable.ExperimentTime, expTimes) & RawDataTable.Location == loc;
    tempRawData = RawDataTable(idx2,:);
    
   tempRawData.Intensity(tempRawData.Channel == 1) = tempRawData.Intensity(tempRawData.Channel == 1) .* redAuToJ(loc);%./ redThresholds(loc);
   tempRawData.Intensity(tempRawData.Channel == 2) = tempRawData.Intensity(tempRawData.Channel == 2) .* greenAuToJ(loc);%./  greenThresholds(loc);
    
    tempRawData = sortrows(tempRawData, 'Intensity');

    numTrials = size(tempRawData,1);
    onesVec = ones(numTrials, 1);
    Red680 = tempRawData.Channel == 1 & tempRawData.Color == 1;
    Green680 = tempRawData.Channel == 1 & tempRawData.Color == 2;
    Achrom680 = tempRawData.Channel == 1 & tempRawData.Color == 3;
    NotSeen680 = tempRawData.Channel == 1 & tempRawData.Color == 4;
    
    Red543 = tempRawData.Channel == 2 & tempRawData.Color == 1;
    Green543= tempRawData.Channel == 2 & tempRawData.Color == 2;
    Achrom543 = tempRawData.Channel == 2 & tempRawData.Color == 3;
    NotSeen543 = tempRawData.Channel == 2 & tempRawData.Color == 4;

    Red680CumSum = cumsum(double(Red680));
    Green680CumSum = cumsum(double(Green680));
    Achrom680CumSum = cumsum(double(Achrom680));
    NotSeen680CumSum = cumsum(double(NotSeen680));

    Red543CumSum = cumsum(double(Red543));
    Green543CumSum = cumsum(double(Green543));
    Achrom543CumSum = cumsum(double(Achrom543));
    NotSeen543CumSum = cumsum(double(NotSeen543));

    CumSumNumTrials = Red680CumSum + Green680CumSum + Achrom680CumSum + NotSeen680CumSum + Red543CumSum + Green543CumSum + Achrom543CumSum + NotSeen543CumSum;
    figure; hold on; title(num2str(propL(loc)));
    plot(tempRawData.Intensity(NotSeen680), NotSeen680CumSum(NotSeen680)./ numTrials, 'Marker', 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'k');
plot(tempRawData.Intensity(Achrom680), Achrom680CumSum(Achrom680)./ numTrials, 'Marker', 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
    plot(tempRawData.Intensity(Red680), Red680CumSum(Red680)./ numTrials, 'Marker', 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    plot(tempRawData.Intensity(Green680), Green680CumSum(Green680)./ numTrials, 'Marker', 'o', 'Color', 'r', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g');


plot(tempRawData.Intensity(NotSeen543), NotSeen543CumSum(NotSeen543)./ numTrials, 'Marker', 'o', 'Color', 'g', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'k');
plot(tempRawData.Intensity(Achrom543), Achrom543CumSum(Achrom543)./ numTrials, 'Marker', 'o', 'Color', 'g', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'w');
plot(tempRawData.Intensity(Red543), Red543CumSum(Red543)./numTrials, 'Marker', 'o', 'Color', 'g', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'r');
    plot(tempRawData.Intensity(Green543), Green543CumSum(Green543)./ numTrials, 'Marker', 'o', 'Color', 'g', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');

set(gca, 'xscale', 'log');

    
end