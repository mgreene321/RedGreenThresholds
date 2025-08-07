function plotStaircase_RG
if isvarname('exp_data')
else
    [dataFile, dataPath] = uigetfile('Z:\Local_Share\Imaging_Data\AOVIS\20256');
    load(fullfile(dataPath, dataFile));
end
data = exp_data.data_matrix;

numLocs = max(data(:,1));
numStair = max(data(:,2));
numAOMs = max(data(:,3));

m = floor(sqrt(numLocs));
n = ceil(sqrt(numLocs));
colors = {'r', 'g'};
for stair = 1:numStair
    figure; hold on
    
    for loc = 1:numLocs
        subplot(m,n,loc); hold on; yticks(0:0.2:1); set(gca, 'yscale', 'log');
        xlabel('Trial')
        ylabel('Stimulus intensity (au)')
        ylim([0 1]);
        
        for aom = 1:numAOMs
            offset = aom-1;
            trialIdx = data(:,1) == loc &...
                data(:,2) == stair &...
                data(:,3) == aom;
            
            nTrials = sum(trialIdx);
            stimVals = data(trialIdx, 4);
            responses = data(trialIdx,end);
            plot(1 + offset:nTrials + offset, stimVals,'Marker', 'o', 'MarkerEdgeColor', colors{aom}, 'Color', colors{aom}, 'LineWidth', 2)
            
            for t = 1:nTrials
                
                if responses(t) == 1
                    scatter(t + offset, stimVals(t), 'ko', 'MarkerFaceColor', 'r')
                elseif responses(t) == 2
                    scatter(t + offset, stimVals(t), 'ko', 'MarkerFaceColor', 'g')
                elseif responses(t) == 3
                    scatter(t + offset, stimVals(t), 'ko', 'MarkerFaceColor', 'w')
                elseif responses(t) == 4
                    scatter(t + offset, stimVals(t), 'ko', 'MarkerFaceColor', 'k')
                end
                
            end
            
            
        end
    end
    
end
end



