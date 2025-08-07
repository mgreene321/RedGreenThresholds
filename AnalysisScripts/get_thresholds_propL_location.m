function ThresholdsPropLstimLocs = get_thresholds_propL_location(subjectId, expDate)

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
subjectFolder = fullfile(RedGreenThresholdsPath, subjectId);
saveFolder = fullfile(subjectFolder, expDate);

analysisDir = dir(fullfile(saveFolder, 'analysisOutput*'));

for i = 1:size(analysisDir,1)
    load(fullfile(analysisDir(i).folder, analysisDir(i).name))
    thresholds{i} = analysisOutput.thresholds;
    
    redThresholdPowerWatts{i}
    
    for j = 1:size(analysisOutput.expTime,1)
        load(fullfile(saveFolder, analysisOutput.expTime{j}, 'LMS.mat'))
        %load(fullfile(saveFolder, analysisOutput.expTime{j}, 'stimLocs.mat'))
        load(fullfile(saveFolder, analysisOutput.expTime{j}, 'deliveryStatistics.mat'), 'meanDeliveryLocs');
        tempL(:,j) = LMS.L;
        tempM(:,j) = LMS.M;
        tempS(:,j) = LMS.S;
        
        
        %tempPropL(:,j) = LMS.L./(LMS.L + LMS.M + LMS.S);
        %tempStimLocs(:,:,j) = StimLocs;
        tempDeliveryLocs(:,:,j) = meanDeliveryLocs;
    end
    %tempPropL = mean(tempPropL, 2);
   % propL{i} = tempPropL;
   
   tempL = mean(tempL,2); tempM = mean(tempM,2); tempS = mean(tempS,2);
   L{i} = tempL; M{i} = tempM; S{i} = tempS;
   
   % stimLocs{i} = mean(tempStimLocs,3); 
   
   deliveryLocs{i} = mean(tempDeliveryLocs,3);
    %tempPropL = [];
    tempL = []; tempM = []; tempS = [];
    
   tempDeliveryLocs = [];
    
    
 
end

thresholds = vertcat(thresholds{:});
    %propL = vertcat(propL{:});
    L = vertcat(L{:});
    M = vertcat(M{:});
    S = vertcat(S{:});
    
    %stimLocs = vertcat(stimLocs{:});
    deliveryLocs = vertcat(deliveryLocs{:});
    
    ThresholdsPropLstimLocs = [thresholds L M S deliveryLocs];

end