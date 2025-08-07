% Make figure and define some LineSpec

figure; hold on
xlabel('Proportion L');
ylabel('S_{680}/S_{543}');
marker = 'o';
markerEdgeColor = 'k';
markerFaceColors = {'k', 'w'}; % 10001R data shown as filled markers, 20217R as open
lineWidth = 2;


%% Plot model predictions

%'propLFine' is a vector of 1000 linearly spaced elements ranging from 0 to
%1, and is the domain over which predictions are made from the additive
%model. RG_XXXnm is the ratio of 680 nm to 543 nm sensitivity derived from
%the additive model, assuming a lambda_max of XXX nm. 

load('additiveModel.mat');

% Plot additive model for lambda_max = 555.5 nm
plot(propLFine, RG_555nm, 'LineStyle', ':', 'LineWidth', lineWidth,...
    'Color', 'k');

% Plot additive model for lambda_max = 558.9 nm, the PTB default
plot(propLFine, RG_559nm, 'LineStyle', '-', 'LineWidth', lineWidth,...
    'Color', 'k');

% Plot additive nmodel for lambda_max = 563.4 nm
plot(propLFine, RG_563nm, 'LineStyle', '--', 'LineWidth', lineWidth,...
    'Color', 'k');


%% Plot subjects' data

% 'tbl' is a table with the following columns:
%   'SubjectID'
%   'PropL': optically weighted proportion L, not counting missing cones
%   'minPropL': proportion L assuming all missing cones are M
%   'maxPropL': proportion L assuming all missing cones are L
%   'RG_ratio': 543 nm corneal irradiance (quanta/s/m^2) at threshold
%   divided by 680 nm corneal irradiance at threshold. 

tbl =importdata('RedGreenThresholdsDataForWill.mat');
subjects = unique(tbl.SubjectID);

for s = 1:numel(subjects)
    subjectIdx = strcmpi(tbl.SubjectID, subjects{s});
    propL = tbl.PropL(subjectIdx);
    minPropL = tbl.minPropL(subjectIdx); maxPropL = tbl.maxPropL(subjectIdx);
    RG_ratio = tbl.RG_ratio(subjectIdx);
    
    errorbar(propL, RG_ratio, [], [], abs(propL-minPropL), abs(propL-maxPropL),...
        'LineStyle', 'none', 'Color', 'k', 'LineWidth', lineWidth, 'Marker', 'none')
    
    plot(propL, RG_ratio, 'LineStyle', 'none', 'Marker', marker,...
        'MarkerEdgeColor', markerEdgeColor, 'MarkerFaceColor',...
        markerFaceColors{s}, 'LineWidth', lineWidth);
end