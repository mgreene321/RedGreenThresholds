clear all;
addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\AOMcontrol\mQUESTPlus-master'));

numLocs = 4;
numStair = 1;
numChrom = 2;
trialsPerStair = 20;
nTrials = numLocs * numStair * trialsPerStair * numChrom;
A = allcomb(1:numLocs, 1:numStair, 1:numChrom);
stim_matrix = repelem(A, trialsPerStair,1);

qpPF = @qpPFNormal;

nOutcomes = 2;
stimParamsDomainList = {0:0.001:1}; %for cum norm
psiParamsDomainList  = {0:0.001:1, 0.1, 0.01}; %SD and lapse are fixed 
stopRule = 'nTrials';

tic;
for loc = 1:numLocs
    for strcs = 1:numStair
        for chrom = 1:numChrom
            
            q = qpInitialize('qpPF', qpPF,...
                'nOutcomes', nOutcomes,...
                'stimParamsDomainList', stimParamsDomainList,...
                'psiParamsDomainList', psiParamsDomainList,...
                'stopRule', stopRule);
            
            questData(loc, strcs, chrom) = q;
        end
    end
end

%Set up simulated observer
simulatedPsiParams = cell(numLocs, numChrom);
for i = 1:numLocs
    for j = 1:numChrom
        simulatedPsiParams{i,j} = [0.1 + 0.8*rand, 0.1, 0.01];
    end
end

simulatedObserverFun = @(x, loc, chrom) qpSimulatedObserver(x, @qpPFNormal, simulatedPsiParams{loc,chrom});

for tt = 1:nTrials
    loc = stim_matrix(tt,1);
    stair = stim_matrix(tt,2);
    aom = stim_matrix(tt,3);
    
    qpIntensity = qpQuery(questData(loc, stair, aom));
    test_intensity = min([max([qpIntensity 0]) 1]);
    
    outcome = simulatedObserverFun(test_intensity, loc, aom);
    
    questData(loc, stair, aom) = qpUpdate(questData(loc, stair, aom), test_intensity, outcome);
end

toc;
%% Plotting
figure; hold on;
lineStyles = {'-', '--'};
lineColors = {'r', 'g'};
markerFills = {'r', 'g', 'w', 'k'};
markerShapes = {'o', 's'};

for loc = 1:numLocs
    subplot(sqrt(numLocs), sqrt(numLocs), loc); hold on
    xlabel('Intensity (au)');
    ylabel('Proportion seen');
    title(['Location ' num2str(loc)]);
    for stair = 1:numStair
        for color = 1:numChrom
            q = questData(loc, stair, color);
            psiParamsIndex = qpListMaxArg(q.posterior);
            psiParamsQuest = q.psiParamsDomain(psiParamsIndex,:);
            psiParamsFit = qpFit(q.trialData, q.qpPF, psiParamsQuest, q.nOutcomes,...
                'lowerBounds', [0 0.01 0.01], 'upperBounds', [1 1 0.01]);
            stimCounts = qpCounts(qpData(q.trialData), q.nOutcomes);
            stim = [stimCounts.stim];
            stimFine = linspace(0,1,100)';
            plotProportionsFit = qpPFNormal(stimFine, psiParamsFit);
            for cc = 1:length(stimCounts)
                nTrials(cc) = sum(stimCounts(cc).outcomeCounts);
                pCorrect(cc) = stimCounts(cc).outcomeCounts(2)/nTrials(cc);
            end
            for cc = 1:length(stimCounts)
                h = scatter(stim(cc),pCorrect(cc),100,'Marker', markerShapes{color},'MarkerEdgeColor','none','MarkerFaceColor',lineColors{color},...
                    'MarkerFaceAlpha',nTrials(cc)/max(nTrials),'MarkerEdgeAlpha',nTrials(cc)/max(nTrials));
            end
            plot(stimFine,plotProportionsFit(:,2),'-','Color',lineColors{color},'LineWidth',3);
%             xlabel('Stimulus Value');
%             ylabel('Proportion Correct');
            
            
        end
    end
end

