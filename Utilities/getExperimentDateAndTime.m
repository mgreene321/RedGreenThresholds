function [dataPaths, expDate, expTime, subjectId] = getExperimentDateAndTime

tempFolder = 'Z:\Local_Share\Imaging_Data\AOVIS';
selectFiles = 1;
i = 1;
while selectFiles
     tempFolder = uigetdir(tempFolder);
    selectFiles = ischar(tempFolder);
    if ~selectFiles; break; end
    dataPaths{i,:} = tempFolder;
    i = i + 1;
end

expDate = regexp(dataPaths , '\d{1,2}_\d{1,2}_\d{4}', 'Match'); % gets mm_dd_yyyy
expTime = regexp(dataPaths, '\d{1,2}_\d{1,2}_\d{1,2}\>', 'Match');
subjectId = regexp(dataPaths,  '\d{5}[LR]|\d{5}.[LR]', 'Match');

expDate = [expDate{:}]';
expTime = [expTime{:}]';

subjectId = cellfun(@(x) erase(x, filesep), subjectId, 'UniformOutput', false); % get rid of \
subjectId = [subjectId{:}]';

if length(unique(expDate)) > 1
    warning('Data from multiple days selected')
end

if length(unique(subjectId)) > 1
    warning('Data from multiple subjects selected')
end

end