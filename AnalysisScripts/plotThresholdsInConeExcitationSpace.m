compareThresholdsToPrediction;
fontSize = 12;
for s = 1:length(subjects)
    subjectId  = subjects{s};
    figure; hold on, set(gca, 'FontSize', fontSize)
    title(subjectId);
    xlabel('L-cone excitations')
    ylabel('M-cone excitations')
    
    subjectId = subjects{s};
    DT = analyzedDataTable(strcmpi(analyzedDataTable.SubjectID, subjectId),:);
    
    propL = DT.L./(DT.L + DT.M);
    
   L_red = DT.L .* effectiveRedQuantaLMS{s}(:,1);
   L_green = DT.L.* effectiveGreenQuantaLMS{s}(:,1);
   
   M_red = DT.M .* effectiveRedQuantaLMS{s}(:,2);
   M_green = DT.M .* effectiveGreenQuantaLMS{s}(:,1);
    
  numSeen543 = DT.NumRed543 + DT.NumGreen543 + DT.NumAchrom543;
  numSeen680 = DT.NumRed680 + DT.NumGreen680 + DT.NumAchrom680;
    
    for i = 1:size(DT,1)
        
%         if propL(i) <= 0.2
%             markerSize = 4;
%         elseif propL(i) <= 0.4
%             markerSize = 5;
%         elseif propL(i) <= 0.6
%             markerSize = 6;
%         elseif propL(i) <= 0.8
%             markerSize = 7;
%         elseif propL(i) <= 1
%             markerSize = 8;
%         end
        
        plot(L_red(i), M_red(i),'Marker', 's', 'MarkerFaceColor', [1 - DT.NumGreen680(i)./numSeen680(i), 1 - DT.NumRed680(i)./numSeen680(i), DT.NumAchrom680(i)./numSeen680(i)], 'MarkerEdgeColor', [propL(i), 1 - propL(i), 0], 'MarkerSize', 8, 'LineWidth', 2)% 'MarkerSize', max(10.* propL(i), 2))
        
        plot(L_green(i), M_green(i),'Marker', 'o', 'MarkerFaceColor', [1 - DT.NumGreen543(i)./numSeen543(i), 1 - DT.NumRed543(i)./numSeen543(i) DT.NumAchrom543(i)./numSeen543(i)],'MarkerEdgeColor', [propL(i),1 - propL(i),0], 'MarkerSize', 8, 'LineWidth', 2) % 'MarkerSize', max(10.* propL(i), 2))
        
        
    end
    
    
    
    
end