root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
%Load analyzed data table
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);
% load raw data
RDT = importdata(fullfile(RedGreenThresholdsPath, 'RawDataTable.mat'));
RDT = RDT(RDT.StimDurFrames == 3,:);

subjects = unique(ADT.SubjectID);
% for each subject

% for each channel (wavelength)

% get probability of hue response for each wavelength
lineWidth = 2;
figure; hold on

xpositions = [1 2; 3.5 4.5];
for s = 1:numel(subjects)
    for ch = 1:2
        for color = 1:3
            PcolorGivenCh{s}{ch}{color} = sum(strcmpi(RDT.SubjectID, subjects{s}) & RDT.Channel == ch & RDT.Color == color)./...
                sum(strcmpi(RDT.SubjectID, subjects{s}) & RDT.Channel == ch & RDT.YesNo == 1);
        end

    end

X = categorical({'543 nm', '680 nm'});
Y(1,:) = horzcat(PcolorGivenCh{s}{2}{:});
Y(2,:) = horzcat(PcolorGivenCh{s}{1}{:});
Y = fliplr(Y);
Y(:,2:end) = fliplr(Y(:,2:end));

X = xpositions(s,:);
b{s} = bar(X, Y, 'stacked', 'LineWidth',2);
set(gca, 'ColorOrder', [1 1 1; 1 0 0; 0 1 0], 'linewidth', lineWidth);
set(gca, 'TickDir', 'out')

end

xticks(sort(xpositions(:)));
xticklabels({'543 nm', '680 nm', '543 nm', '680 nm'})
set(gca, 'fontsize', 16)
