root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);

fontSize = 16;
figure; hold on; set(gca, 'Position', [0.175 0.175 0.8 0.8], 'FontSize', fontSize, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 3);

xlabel('Heterogeneity Index')
ylabel('Proportion color response');
markerFaceColors = {'k', 'w'};

subjects = unique(ADT.SubjectID);
propL = ADT.L./(ADT.L + ADT.M);
minPropL = ADT.L./(ADT.L + ADT.M + ADT.X);
maxPropL = (ADT.L + ADT.X)./(ADT.L + ADT.M + ADT.X);


heteroIdx = 1 - 2.*abs(0.5 - propL);
heteroIdxMinPropL = 1 - 2.*abs(0.5-minPropL);
heteroIdxMaxPropL = 1 - 2.*abs(0.5-maxPropL);

heteroIdxXneg = nan(size(heteroIdx));
heteroIdxXpos = nan(size(heteroIdx));
for i = 1:size(heteroIdx,1)
    if heteroIdxMinPropL(i) < heteroIdx(i) & heteroIdx(i) < heteroIdxMaxPropL(i)
        heteroIdxXneg(i) = abs(heteroIdx(i) - heteroIdxMinPropL(i));
        heteroIdxXpos(i) = abs(heteroIdx(i) - heteroIdxMaxPropL(i));
    elseif heteroIdxMaxPropL(i) < heteroIdx(i) & heteroIdx(i) < heteroIdxMinPropL(i)
         heteroIdxXneg(i) = abs(heteroIdx(i) - heteroIdxMaxPropL(i));
        heteroIdxXpos(i) = abs(heteroIdx(i) - heteroIdxMinPropL(i));
        
    elseif heteroIdxMinPropL(i) < heteroIdx(i) & heteroIdxMaxPropL(i) < heteroIdx(i)
        heteroIdxXneg(i) = abs(heteroIdx(i) - min([heteroIdxMinPropL(i) heteroIdxMaxPropL(i)]));
    elseif heteroIdxMinPropL(i) > heteroIdx(i) & heteroIdxMaxPropL(i) > heteroIdx(i)
        heteroIdxXpos(i) = abs(heteroIdx(i) - max([heteroIdxMinPropL(i) heteroIdxMaxPropL(i)]));
    end
    
end



% Color responses out of seen responses

colorOutOfSeen = (ADT.NumRed680 + ADT.NumGreen680 + ADT.NumRed543 + ADT.NumGreen543)...
    ./(ADT.NumRed680 + ADT.NumGreen680 + ADT.NumRed543 + ADT.NumGreen543 + ADT.NumAchrom543 + ADT.NumAchrom680);

for s = 1:length(subjects)
subjectId = subjects{s};
subjectIdx = strcmpi(ADT.SubjectID, subjectId); 

errorbar(heterogeneityIdx(subjectIdx), colorOutOfSeen(subjectIdx), [], [], heteroIdxXneg(subjectIdx), heteroIdxXpos(subjectIdx), ...
        'LineStyle', 'none', 'LineWidth', 2, 'Marker', 'none', 'Color', 'k');


plot(heterogeneityIdx(subjectIdx), colorOutOfSeen(subjectIdx),...
    'LineStyle', 'none', 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', markerFaceColors{s}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'DisplayName', subjects{s})

end


% bin data for each subject
binEdges = 0:(1/7):1;
midPoints = binEdges(1:end-1) + diff(binEdges)/2;


for subj = 1:numel(subjects)
    subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
    colorOutOfSeenSubj = colorOutOfSeen(subjIdx);
    [N{subj}, ~, bin{subj}] = histcounts(heterogeneityIdx(subjIdx), binEdges);
    
    for b = 1:numel(midPoints)
        
        colorOutOfSeenBinned{subj}(b) = mean(colorOutOfSeenSubj(bin{subj} == b));
        
    end
    
end