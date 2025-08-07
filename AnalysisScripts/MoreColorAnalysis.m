clear;
root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);


lineWidth = 2;
markerSize = 8;
fontSize = 12;

A = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
subjectId = '20217R';
A = A(strcmpi(A.SubjectID, subjectId),:);

propL = A.L./(A.L + A.M);

red = A.NumRed680 + A.NumRed543;
green = A.NumGreen680 + A.NumGreen543;
achrom = A.NumAchrom680 + A.NumAchrom543;
notSeen = A.NumNotSeen680 + A.NumNotSeen543;

seen680 = A.NumRed680 + A.NumGreen680 + A.NumAchrom680;
notSeen680 = A.NumNotSeen680;

seen543 = A.NumRed543 + A.NumGreen543 + A.NumAchrom543;
notSeen543 = A.NumNotSeen543;
 
seen = seen680 + seen543;

%% Proportion of non-achromatic seen responses vs. proportion L

propChromOutOfSeen = (red + green)./(red + green + achrom);

propRedOutOfSeen = red./(red + green + achrom);
propGreenOutOfSeen = green./(red+green+achrom);


% bin data
edges = 0:0.2:1;
midPoints = edges(1:end-1) + diff(edges)/2;

[N,~, bin] = histcounts(propL, edges);

propChromOutOfSeenBinned = zeros(size(midPoints));
propChromOutOfSeenSeenError = zeros(size(midPoints));

propRedOutOfSeenBinned = zeros(size(midPoints));
propGreenOutOfSeenBinned = zeros(size(midPoints));



for i = 1:length(midPoints)
    propChromOutOfSeenBinned(i) = mean(propChromOutOfSeen(bin==i), 'omitnan');
    propChromOutOfSeenSeenError(i) = std(propChromOutOfSeen(bin==i), 'omitnan');
    
    propRedOutOfSeenBinned(i) =  mean(propRedOutOfSeen(bin==i), 'omitnan');
    propRedOutOfSeenError(i) = std(propRedOutOfSeen(bin==i), 'omitnan');
    
    propGreenOutOfSeenBinned(i) = mean(propGreenOutOfSeen(bin==i), 'omitnan');
    propGreenOutOfSeenError(i) = std(propGreenOutOfSeen(bin==i), 'omitnan');
    
end

%figure(1); hold on;errorbar(midPoints,propRedOutOf680Binned, propRedOutOf680Error, 'r');errorbar(midPoints,propRedOutOf680Binned, propRedOutOf680Error, 'r');
% figure; hold on
% 
% plot(midPoints, propChromOutOfSeenBinned, 'ko-')
% errorbar(midPoints, propChromOutOfSeenBinned, propChromOutOfSeenSeenError, 'k', 'LineWidth', lineWidth)
% 
% plot(midPoints, propRedOutOfSeenBinned, 'ro-', 'LineWidth', lineWidth);
% errorbar(midPoints, propRedOutOfSeenBinned, propRedOutOfSeenError, 'r', 'LineWidth', lineWidth);
% 
% plot(midPoints, propGreenOutOfSeenBinned, 'go-', 'LineWidth', lineWidth);
% errorbar(midPoints, propGreenOutOfSeenBinned, propGreenOutOfSeenError, 'g', 'LineWidth', lineWidth);

%% Breakdown of responses to 680 or 543 nm vs proportion L

% 680 nm

propRedOutOf680 = A.NumRed680./seen680;
propGreenOutOf680 = A.NumGreen680./seen680;
propAchromOutOf680 = A.NumAchrom680./seen680;


% 543

propRedOutOf543 = A.NumRed543./seen543;
propGreenOutOf543 = A.NumGreen543./seen543;
propAchromOutOf543 = A.NumAchrom543./seen543;

for i = 1:length(midPoints)
    % 680
    propRedOutOf680Binned(i) = mean(propRedOutOf680(bin==i), 'omitnan');
    propRedOutOf680Error(i) = std(propRedOutOf680(bin==i), 'omitnan');
    
    propGreenOutOf680Binned(i) = mean(propGreenOutOf680(bin==i), 'omitnan');
    propGreenOutOf680Error(i) = std(propGreenOutOf680(bin==i), 'omitnan');
    
%     propAchromOutOf680Binned(i) = mean(propAchromOutOf680(bin==i), 'omitnan');
%     propAchromOutOf680Error(i) = std(propAchromOutOf680(bin==i), 'omitnan');
    
    % 543
    propRedOutOf543Binned(i) = mean(propRedOutOf543(bin==i), 'omitnan');
    propRedOutOf543Error(i) = std(propRedOutOf543(bin==i), 'omitnan');
    
    propGreenOutOf543Binned(i) = mean(propGreenOutOf543(bin==i), 'omitnan');
    propGreenOutOf543Error(i) = std(propGreenOutOf543(bin==i), 'omitnan');
    
%     propAchromOutOf543Binned(i) = mean(propAchromOutOf543(bin==i), 'omitnan');
%     propAchromOutOf543Error(i) = std(propAchromOutOf543(bin==i), 'omitnan');

end

figure; hold on; title(['Red (680 nm) Stimuli, ' subjectId]), set(gca, 'FontSize', fontSize);
xlim([0 1]); xlabel('Proportion L')
ylim([0 1]); ylabel('Response Proportion')

for i = 1:numel(midPoints)
plot(midPoints(i), propRedOutOf680Binned(i), 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', N(i)+1, 'LineWidth', lineWidth)
end
errorbar(midPoints,propRedOutOf680Binned, propRedOutOf680Error, 'Color', 'r', 'LineWidth', lineWidth, 'MarkerSize', markerSize);

for i = 1:numel(midPoints)
plot(midPoints(i), propGreenOutOf680Binned(i), 'Marker', 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', N(i)+1, 'LineWidth', lineWidth)
end
errorbar(midPoints,propGreenOutOf680Binned, propGreenOutOf680Error, 'Color', 'g', 'LineWidth', lineWidth, 'MarkerSize',  markerSize);
%plot(midPoints, propAchromOutOf680Binned, 'ko-')
%plot(midPoints, propNotSeenOutOf680Binned, 'ko--')


figure; hold on; title(['Green (543 nm) Stimuli, ' subjectId]), set(gca, 'FontSize', fontSize);
xlim([0 1]); xlabel('Proportion L')
ylim([0 1]); ylabel('Response Proportion')

for i = 1:numel(midPoints)
plot(midPoints(i), propRedOutOf543Binned(i),'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r',  'MarkerSize', N(i)+1, 'LineWidth', lineWidth)
end
errorbar(midPoints,propRedOutOf543Binned, propRedOutOf543Error, 'Color', 'r', 'LineWidth', lineWidth, 'MarkerSize',  markerSize);

for i = 1:numel(midPoints)
plot(midPoints(i), propGreenOutOf543Binned(i),'Marker', 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g',  'MarkerSize', N(i)+1, 'LineWidth', lineWidth)
end
errorbar(midPoints,propGreenOutOf543Binned, propGreenOutOf543Error, 'Color', 'g', 'LineWidth', lineWidth, 'MarkerSize',  markerSize);



%plot(midPoints, propAchromOutOf543Binned, 'ko-')
%plot(midPoints, propNotSeenOutOf543Binned, 'ko--')








