function RedGreenThresholds_Max


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function runs an experiment for measuring detection thresholds.
%
%  02/06/19 JEV - Wrote it.
%  06/12/21 MG - Edited it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
rng('shuffle');
% set some variables to global; most of these are first modified
% by AOMcontrol.m
global SYSPARAMS StimParams VideoParams;

addpath(genpath('C:\Programs\AOMcontrol\Experiments\mQUESTPlus-master'));
%addpath(genpath('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds'));

%% Hardcoded parameters
numColor = 2; %number of stimulus chromaticities is 2 (just red or green)
videodur = 1;

%% Set keys

valid_keys = {'space', 'shift', 'rightarrow', 'leftarrow', 'uparrow', 'downarrow', 'escape'};
[keys.startTrial, keys.repeatTrial,...
    keys.stimRed, keys.stimGreen,...
    keys.stimAchromatic, keys.stimNotSeen,...
    keys.quitKey] = deal(valid_keys{:});


%% Load UI
[sessionParams, expParams] = loadGUI;
numLocs = expParams.gridWidth * expParams.gridHeight;
numDiams = numel(expParams.stimDiam);
numDurs = numel(expParams.stimDur);


%See set_VideoParams_PsyfileName.m
filename = fullfile(VideoParams.rootfolder,...
    ['RedGreenThresholds_',sessionParams.initials, '_',strrep(strrep(strrep(datestr(now), ...
    '-',''), ' ','x'),':',''),'.mat']); %DDMMMYYYY

sessionParams.filename = filename;

%% Compute locations

[X Y] = meshgrid(0:expParams.gridHeight - 1, 0:expParams.gridWidth - 1);
stimOff = max(expParams.stimDiam) .* [X(:) Y(:)];
A = allcomb(1:numLocs, 1:expParams.numStair, 1:numColor, 1:numDiams, 1:numDurs);

%% Stimulus matrix
%Col 1: Location index
%Col 2: Staircase number
%Col 3: Color/aom index

%nTrials = numLocs * expParams.numStair * expParams.trialsPerStair * numColor;

%columns: loc stair color diam dur
stim_matrix = repelem(A, expParams.trialsPerStair,1);
%stim_matrix(:,end+1) = rand(size(stim_matrix,1),1);
%stim_matrix = sortrows(stim_matrix, size(stim_matrix,2));


% data matrix is stim mat with 3 columns of nan for stim intensity, seen
% flag and color percept respectively

%Add in refresh trials

refreshTrialProportion = 0.2;
refreshTrialsPerStair = round(expParams.trialsPerStair./((1/refreshTrialProportion)-1));
refreshTrial_matrix = repelem(A, refreshTrialsPerStair, 1);
refreshTrial_matrix(:,2) = 0; %set staircase number to 0 to mark refresh trials

stim_matrix = [stim_matrix; refreshTrial_matrix];
stim_matrix(:,end+1) = rand(size(stim_matrix,1),1);
stim_matrix = sortrows(stim_matrix, size(stim_matrix,2));

%columns: columns: loc stair color diam dur rand
exp_data.stim_matrix = stim_matrix;
exp_data.data_matrix = [stim_matrix(:,1:end-1) nan(size(stim_matrix,1),3)];

exp_data.sessionParams = sessionParams;
exp_data.expParams = expParams;
save(filename,'exp_data');

nTrials = size(stim_matrix,1);

% figure out the trial number of the first trial of each qp staircase
firstTrials = [];

for loc = 1:numLocs
    for diam = 1:numDiams
        for dur = 1:numDurs
            for stair = 1:expParams.numStair
                for color = 1:numColor
                    
                    trialsOfInterest = find(stim_matrix(:,1) == loc & stim_matrix(:,2) == stair & stim_matrix(:,3) == color...
                        & stim_matrix(:,4) == diam & stim_matrix(:,5) == dur);
                    
                    firstTrials = [firstTrials trialsOfInterest(1)];
                    
                end
            end
        end
    end
end


%After loading custom reference

%% Netcomm commands

% Set and apply transverse chromatic offsets
setTCA = ['UpdateTCA#' num2str(expParams.redTCA_X) '#'...
    num2str(expParams.redTCA_Y) '#'...
    num2str(expParams.greenTCA_X) '#'...
    num2str(expParams.greenTCA_Y) '#'...
    '0#0#'];
applyTCA = 'ApplyTCA#1#';

% Set duty cycle to 1:0
setFlash = 'Flash#30#0#';

% Set gain and turn on proper digital marks
gain = expParams.gain;
setGain = ['Gain#' num2str(gain) '#'];
digitalMarks = ['DigitalMarks' repmat(['#' num2str(gain)], 1,3)...
    repmat(['#' num2str(~gain)], 1,3) '#'];

% Set stimulus location
setGridOrigin = ['LocUpdateAbs#' num2str(expParams.targetX) '#'...
    num2str(expParams.targetY), '#'];

% Commmunicate subject ID
setPrefix = ['VP#' sessionParams.initials '#'];

% Deliver commands
commands = {digitalMarks, setTCA, applyTCA, setFlash, setGain, setGridOrigin, setPrefix};

if SYSPARAMS.realsystem
    for i = 1:length(commands)
        netcomm('write', SYSPARAMS.netcommobj, int8(commands{i})); pause(0.2);
    end
else
end
pause(0.2);

Speak('Apply TCA.');

%% Initialize experiment
% if cone_select.main_gui has been successful, initialize exp. Do this here
% so that if video was not good, then will not update parameters in ICANDI
global aom_fig_handle
[hAomControl, aom_fig_handle] = exp.initialize(sessionParams);

%% Set up Mov structure

for dur = 1:numDurs
    [Mov, aom_seq_generic(dur,:)] = initializeMov(expParams.stimDur(dur), gain);
end
StimParams.fext = 'buf';

%% QUEST+ setup
tMin = 0;
tMax = 1;
%QUEST options
%qpPF = @qpPFWeibull;
qpPF = @qpPFNormal;

nOutcomes = 2;
stimParamsDomainList = {0:0.001:1}; %for cum norm
psiParamsDomainList  = {0:0.001:1,0.05, 0.01}; %mean, SD are fixed
stopRule = 'nTrials';

q = qpInitialize('qpPF', qpPF,...
    'nOutcomes', nOutcomes,...
    'stimParamsDomainList', stimParamsDomainList,...
    'psiParamsDomainList', psiParamsDomainList,...
    'stopRule', stopRule);

for loc = 1:numLocs
    for diam = 1:numDiams
        for dur = 1:numDurs
            for strcs = 1:expParams.numStair
                for aom = 1:numColor
                    q.stimOffsets = stimOff(loc,:);
                    q.crossLoc = [expParams.targetX expParams.targetY];
                    q.aom = aom;
                    q.stimDiam = expParams.stimDiam(diam);
                    q.stimDur = expParams.stimDur(dur);
                    questData(loc, strcs, aom, diam, dur) = q;
                end
            end
        end
    end
end

%% Create default stimulus
%Note that createStimulus also creates an IR decrement stimulus called
%frame2.bmp
createStim(expParams.stimDiam(1), 'square', 1, 'buf');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Experiment loop                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global runExperiment
runExperiment = 1;
state = 'a';
trial = 1;
resp_flag = 0;
respConfirmed = 0;
set(aom_fig_handle.aom_main_figure, 'KeyPressFcn','uiresume');

Speak('Begin experiment');

while(runExperiment == 1)
    %---------------- State a: quit or start ----------------------------------
    if state == 'a'                                                        % Terminate experiment if maximum trials reached
        if trial > nTrials
            Speak('Experiment complete!')
            save(filename, 'questData', '-append');
            runExperiment = 0;
            uiresume;
            TerminateExp;
            message = ['Off - Experiment Aborted - Trial ' num2str(trial) ...
                ' of ' num2str(nTrials)];
            set(aom_fig_handle.aom1_state, 'String', message);
            resp = [];
        elseif (trial > 1 && trial <= nTrials)                             % After first trial, before final trial, go to state b
            respConfirmed = 0;
            resp = [];
            state = 'b';
        else                                                               % First trial: get response
            uiwait;
            resp = get(aom_fig_handle.aom_main_figure,'CurrentKey');
        end
        
        if strcmpi(resp, keys.quitKey) && state == 'a'                     % Quit on first trial
            % Make sure user really wants to quit
            
            answer = questdlg('Do you really want to quit?', 'Quit confirmation', 'Yes', 'No', 'No');
            
            if strcmpi(answer, 'Yes')
                save(filename, 'questData', '-append');
                abortExperiment(trial, nTrials); break;
                %TerminateExp;
                
            else
                resp = [];
            end
        elseif strcmpi(resp, keys.startTrial)                              % Start first trial, go to state b
            state = 'b';
            
        end
    end
    
    %---------------- State b: determine and play stimulus  -------------------
    if state == 'b'
        
        test_loc = stim_matrix(trial,1); %gets location number
        active_staircase = stim_matrix(trial,2); %gets strcs number
        active_aom = stim_matrix(trial,3); %gets aom number
        diamIdx = stim_matrix(trial, 4); test_diam = expParams.stimDiam(diamIdx);
        durIdx = stim_matrix(trial, 5); test_dur = expParams.stimDur(durIdx);
        
        if active_aom == 1
            aom_color = 'red';
        elseif active_aom == 2
            aom_color = 'green';
        end
        
        if active_staircase == 0 %refresh trial
            test_intensity = round(rand,3);
        else
            q = questData(test_loc, active_staircase, active_aom, diamIdx, durIdx);
            qpIntensity = qpQuery(q);
            
            
            % Modify test intensity according to condition
            if ismember(trial, firstTrials)
                if mod(active_staircase,2) == 1 %odd numbered staircase
                    test_intensity = 0.25;
                elseif mod(active_staircase,2) == 0 % even numbered staircase
                    test_intensity = 0.75;
                end
            else
                test_intensity = min([max([qpIntensity tMin]) tMax]);
            end
        end
        
        
        printTrialInfo(test_intensity, test_loc, active_staircase, aom_color, test_diam, test_dur);
        %disp(qpIntensity);
        
        % play sound to indicate start of stimulus
        beep;
        % change the message displayed in status bar
        message = ['Running Experiment - Trial ' num2str(trial) ...
            ' of ' num2str(nTrials)];
        set(aom_fig_handle.aom1_state, 'String', message);
        Mov.msg = message;
        Mov.seq = '';
        
        % create new stimuli for each trial with new intensity
        createStim(expParams.stimDiam(diamIdx),'square', ...
            test_intensity, 'buf');
        
        if SYSPARAMS.realsystem == 1 %0 on my PC/mac, 1 irl
            StimParams.sframe = 2; %was 2
            StimParams.eframe = 5;
            %StimParams.eframe = size(a,1); %# of buf files
            Parse_Load_Buffers(0);
        end
        
        % ---- set movie parameters to be played by aom ---- %
        % Select AOM power 100% for most experiments unless set
        % otherwise with intensity variable at top of file.
        Mov.aom1pow(:) = 1; % Red (680 nm)
        Mov.aom2pow(:) = 1; % Green (543 nm)
        Mov.aom0pow(:) = 1; % IR (aka imaging)
        
        if strcmpi(aom_color, 'red')
            Mov.aom1offx = stimOff(test_loc,1).*ones(size(aom_seq_generic(durIdx,:)));
            Mov.aom1offy = stimOff(test_loc,2).*ones(size(aom_seq_generic(durIdx,:)));
            Mov.aom0seq = 5.*aom_seq_generic(durIdx,:);
            % Update Mov structure on every trial
            Mov.aom2seq(:) = 0; % switch green off
            Mov.aom1seq = 4.*aom_seq_generic(durIdx,:); % switch red on
        elseif strcmpi(aom_color, 'green')
            
            Mov.aom2offx = stimOff(test_loc,1).*ones(size(aom_seq_generic(durIdx,:)));
            Mov.aom2offy = stimOff(test_loc,2).*ones(size(aom_seq_generic(durIdx,:)));
            Mov.aom0seq = 5.*aom_seq_generic(durIdx,:);
            % Update Mov structure on every trial
            Mov.aom1seq(:) = 0; % switch red off
            Mov.aom2seq = 4.*aom_seq_generic(durIdx,:); % switch green on
        end
        
        % update save name of video
        VideoParams.vidname = [sessionParams.initials '_' sprintf('%03d', trial)];
        
        % send the Mov structure to app data
        setappdata(hAomControl, 'Mov', Mov);
        
        % use the Mov structure to play a movie
        PlayMovie; %%pause(0.1);
        
        state = 'c';
    end
    
    %---------------- State c: get and record response, update staircase,------
    
    if state == 'c'
        
        while respConfirmed == 0
            
            if resp_flag == 0
                uiwait;
                resp = get(aom_fig_handle.aom_main_figure,'CurrentKey');
                disp(resp)
                resp_flag = 1;
            end
            
            if resp_flag == 1
                if ~any(strcmpi(resp, valid_keys))
                    e1 = [Mov.msg ' Invalid response'];
                    set(aom_fig_handle.aom1_state, 'String', e1);
                    resp_flag = 0;
                elseif strcmpi(resp, keys.quitKey)
                    
                    % Make sure user really wants to quit
                    
                    answer = questdlg('Do you really want to quit?', 'Quit confirmation', 'Yes', 'No', 'No');
                    
                    if strcmpi(answer, 'Yes')
                        save(filename, 'questData', '-append');
                        abortExperiment(trial, nTrials); break;
                        %TerminateExp;
                        
                    else
                        resp = [];
                        resp_flag = 0;
                    end
                    
                elseif strcmpi(resp, keys.repeatTrial)  || strcmpi(resp, keys.startTrial)
                    PlayMovie;
                    resp_flag = 0;
                elseif strcmpi(resp, keys.stimNotSeen)
                    seen_flag = 0;
                    color_per = 4;
                    resp_flag = 2;
                elseif strcmpi(resp, keys.stimRed)
                    seen_flag = 1;
                    color_per = 1;
                    resp_flag = 2;
                elseif strcmpi(resp, keys.stimGreen)
                    seen_flag  = 1;
                    color_per = 2;
                    resp_flag = 2;
                elseif strcmpi(resp, keys.stimAchromatic)
                    seen_flag = 1;
                    color_per = 3;
                    resp_flag = 2;
                else
                end
            end
            
            if resp_flag == 2
                uiwait;
                resp = get(aom_fig_handle.aom_main_figure,'CurrentKey');
                disp(resp)
                
                if strcmpi(resp, keys.startTrial)
                    %LOG RESP
                    exp_data.data_matrix(trial, end-2:end) = ...
                        [test_intensity seen_flag color_per];
                    
                    % save(filename,'-append','exp_data');
                    
                    %questData(test_loc, active_staircase, active_aom) = q;
                    %SAVE TRIAL
                    
                    save(filename,'-append','exp_data');
                    
                    if active_staircase ~= 0
                        q  = qpUpdate(q, test_intensity, seen_flag+1);
                        questData(test_loc, active_staircase, active_aom, diamIdx, durIdx) = q;
                    else
                    end
                    
                    %Reset loop variables
                    trial = trial + 1;
                    state = 'a';
                    resp_flag = 0;
                    respConfirmed = 1;
                    
                elseif strcmpi(resp, keys.quitKey)
                    
                    % Make sure user really wants to quit
                    
                    answer = questdlg('Do you really want to quit?', 'Quit confirmation', 'Yes', 'No', 'No');
                    
                    if strcmpi(answer, 'Yes')
                        save(filename, 'questData', '-append');
                        abortExperiment(trial, nTrials); break;
                        %TerminateExp;
                        
                    else
                        resp = [];
                        resp_flag = 0;
                    end
                elseif strcmpi(resp, keys.repeatTrial)
                    PlayMovie;
                    resp_flag = 0;
                    
                elseif ~any(strcmpi(resp, valid_keys))
                    e1 = [Mov.msg ' Invalid response'];
                    set(aom_fig_handle.aom1_state, 'String', e1);
                    resp_flag = 0;
                else
                    resp_flag = 1;
                end
            else
            end
            
        end
    end
end
%% Save and organize data

save(filename, 'questData', '-append');

vidDir = dir(fullfile(VideoParams.rootfolder, sessionParams.initials)); %On ICANDI PC
vidDir = vidDir(find(~cellfun(@isdir,{vidDir(:).name})));

for i = 1:size(vidDir,1)
    if vidDir(i).name(end-3) == '.' %As in '.mat' or '.tif'
        vidDir(i).datenum = 0;
    end
end

[A,I] = max([vidDir(:).datenum]);
latest_folder = vidDir(I).name; %mm_dd_yyyy_HH_MM_SS format

mat_folder = fullfile(VideoParams.rootfolder, sessionParams.initials, latest_folder, filesep);
movefile(filename, mat_folder)

%% PLot

data = exp_data.data_matrix;

numLocs = max(data(:,1));
numStair = max(data(:,2));
numAOMs = max(data(:,3));

m = floor(sqrt(numLocs));
n = ceil(sqrt(numLocs));
colors = {'r', 'g'};
for diam = 1:numDiams
    for dur = 1:numDurs
        for stair = 1:numStair
            figure; hold on
            
            for loc = 1:numLocs
                subplot(m,n,loc); hold on;
                
                title(['Location ' num2str(loc) ' Diam=' num2str(expParams.stimDiam(diam)) ' Dur=' num2str(expParams.stimDur(dur))]);
                ylim([0 1]); set(gca, 'yscale', 'log');
                xlabel('Trial')
                ylabel('Stimulus intensity (au)')
                
                
                for aom = 1:numAOMs
                    offset = aom-1;
                    trialIdx = data(:,1) == loc &...
                        data(:,2) == stair &...
                        data(:,3) == aom & data(:,4) == diam & data(:,5) == dur;
                    
                    nTrials = sum(trialIdx);
                    stimVals = data(trialIdx, end-2);
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
            %    sgtitle(['Staircase ' num2str(stair)])
        end
    end
end

% OPlot refresh trials

for diam = 1:numDiams
    for dur = 1:numDurs
        figure; hold on;
        stair = 0;
        for loc = 1:numLocs
            subplot(m,n,loc); hold on; title(['Refresh, Location ' num2str(loc) ' Diam=' num2str(expParams.stimDiam(diam)) ' Dur=' num2str(expParams.stimDur(dur))]);
            ylim([0 1]); set(gca, 'yscale', 'log');
            xlabel('Trial')
            ylabel('Stimulus intensity (au)')
            for aom = 1:numAOMs
                offset = aom-1;
                trialIdx = data(:,1) == loc &...
                    data(:,2) == stair &...
                    data(:,3) == aom & data(:,4) == diam & data(:,5) == dur;
                
                nTrials = sum(trialIdx);
                stimVals = data(trialIdx, end-2);
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
%sgtitle('Refresh trials')

% plot PFS

try
    
    psiParamsFit = fitAndPlotPFs(data);
    disp(psiParamsFit);
    
catch
end


%% Auxiliary functions

function printTrialInfo(test_intensity, test_loc, active_staircase, aom_color, test_diam, test_dur)

fprintf('\nCurrent intensity: %3f',test_intensity);
fprintf('\nCurrent location: %d',test_loc);
fprintf('\nCurrent staircase: %d',active_staircase);
fprintf('\nCurrent AOM: %s', aom_color);
fprintf('\nCurrennt diameter (pix): %3f', test_diam);
fprintf('\nCurrent duration (f): %3f', test_dur);


function abortExperiment(trial, nTrials)
global runExperiment aom_fig_handle
Speak('Experiment aborted');
runExperiment = 0;
uiresume;
TerminateExp;
message = ['Off - Experiment Aborted - Trial ' num2str(trial) ...
    ' of ' num2str(nTrials)];
set(aom_fig_handle.aom1_state, 'String', message);


function [Mov, aom_seq_generic] = initializeMov(stimDur,gain)
global StimParams
fps = 30;
videodur = 1; %seconds
presentdur_frames = stimDur;
startframe = floor(fps*videodur/2) - floor(presentdur_frames/2);
aom_seq_generic = zeros(1, fps*videodur);
aom_seq_generic(startframe:startframe + presentdur_frames - 1) = 1; %frame num = 1

Mov.duration = fps*videodur;
Mov.gainseq = gain.*ones(1,fps*videodur);
Mov.angleseq = zeros(1,fps*videodur);
Mov.frm = 1;
Mov.seq = '';
Mov.stimbeep = 0.*aom_seq_generic;
Mov.dir = StimParams.stimpath;
Mov.suppress = 0;
Mov.pfx = StimParams.fprefix;

% IR params
Mov.aom0seq = zeros(1,fps*videodur);
Mov.aom0pow = zeros(1,fps*videodur);
Mov.aom0locx = zeros(1,fps*videodur);
Mov.aom0locy = zeros(1,fps*videodur);

% Red params
Mov.aom1seq = zeros(1,fps*videodur);
Mov.aom1pow = zeros(1,fps*videodur);
Mov.aom1offx = zeros(1,fps*videodur);
Mov.aom1offy = zeros(1,fps*videodur);

% Green params
Mov.aom2seq = zeros(1,fps*videodur);
Mov.aom2pow = zeros(1,fps*videodur);
Mov.aom2offx = zeros(1,fps*videodur);
Mov.aom2offy = zeros(1,fps*videodur);

function psiParamsFit = fitAndPlotPFs(goodData)

locIdx = unique(goodData(:,1));
colorIdx = unique(goodData(:,3));
numLocs = length(locIdx);
numDiams = length(unique(goodData(:,4)));
numDurs = length(unique(goodData(:,5)));

%figure; hold on;
lineColors = {'r', 'g'};
markerShapes = {'o', 's'};

%fit weibulls
for diam = 1:numDiams
    for dur = 1:numDurs
        figure; hold on;
        for l = 1:numLocs
            loc = locIdx(l);
            subplot(ceil(sqrt(numLocs)),floor(sqrt(numLocs)),loc); hold on
            title(['Location: ' num2str(loc) ' Diam ' num2str(diam) ' Dur ' num2str(dur)]);
            
            %sgtitle(['Diam ' num2str(diam) 'Dur ' num2str(dur)]);
            for c = 1:2
                color = colorIdx(c);
                tempData = goodData(goodData(:,1) == loc & goodData(:,3) == color & goodData(:,4) == diam & goodData(:,5) == dur,:);
                trialData = [];
                for i = 1:size(tempData,1)
                    trialData(i).stim = 20.*log10(tempData(i,end-2));
                    trialData(i).outcome = tempData(i,end-1) + 1;
                end
                % Set up fake quest+
                qpPF = @qpPFWeibull;
                nOutcomes = 2;
                stimParamsDomainList = {20.*log10(0:0.001:1)}; %for cum norm
                psiParamsDomainList  = {20.*log10(0:0.001:1), 3, 0.01, 0.01}; %mean, SD are fixed
                stopRule = 'nTrials';
                q = qpInitialize('qpPF', qpPF,...
                    'nOutcomes', nOutcomes,...
                    'stimParamsDomainList', stimParamsDomainList,...
                    'psiParamsDomainList', psiParamsDomainList,...
                    'stopRule', stopRule, 'noentropy', true);
                for i = 1:size(tempData,1)
                    q = qpUpdate(q, trialData(i).stim, trialData(i).outcome);
                end
                psiParamsIndex = qpListMaxArg(q.posterior);
                psiParamsQuest = q.psiParamsDomain(psiParamsIndex,:);
                psiParamsFit{l,c} = qpFit(q.trialData, q.qpPF, psiParamsQuest, q.nOutcomes,...
                    'lowerBounds', [20.*log10(0) 0.5 0.01 0.01], 'upperBounds', [20.*log10(1) 20 0.01 0.01]);
                
                % struct with two fields: stim (stimulus intensity) and
                % outcomeCounts (first column gives number of "not seen" responses,
                % second column gives number of "seen" responses, for correspodning
                % stimulus intensity)
                stimCounts = qpCounts(qpData(q.trialData), q.nOutcomes);
                stim = [stimCounts.stim];
                stimFine = 20.*log10(linspace(0,1,1000))';
                plotProportionsFit = qpPFWeibull(stimFine, psiParamsFit{l,c});
                % for each stimulus intensity
                for cc = 1:length(stimCounts)
                    nTrials(cc) = sum(stimCounts(cc).outcomeCounts);
                    pCorrect(cc) = stimCounts(cc).outcomeCounts(2)/nTrials(cc);
                end
                for cc = 1:length(stimCounts)
                    h = scatter(10.^(stim(cc)./20),pCorrect(cc),100,'Marker', markerShapes{c},'MarkerEdgeColor','none','MarkerFaceColor',lineColors{c},...
                        'MarkerFaceAlpha',nTrials(cc)/max(nTrials),'MarkerEdgeAlpha',nTrials(cc)/max(nTrials));
                end
                plot(10.^(stimFine./20),plotProportionsFit(:,2),'-','Color',lineColors{c},'LineWidth',3);
                %             xlabel('Stimulus Value');
                %             ylabel('Proportion Correct');
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                GUI                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sessionParams, expParams] = loadGUI

paramsPath = 'C:\Programs\AOVIS_exp\+Max_exp\RedGreenThresholds';
%Which subject
subjectId = input('Subject ID: ', 's');
fprintf('\n');

%Make subject folder for storing GUI input, if it doesn't already exist
subjectPath = fullfile(paramsPath, subjectId);
if ~isdir(subjectPath)
    mkdir(subjectPath)
else
end

subjectDir = dir(fullfile(subjectPath, '*.mat'));

%Make latest inputs the default
if ~isempty(subjectDir)
    mostRecentIdx = [subjectDir(:).datenum] == max([subjectDir(:).datenum]);
    load(fullfile(subjectPath,subjectDir(mostRecentIdx).name));
    expParamsDef = cellfun(@num2str, struct2cell(expParams), 'UniformOutput', false);
    sessionParamsDef = cellfun(@num2str, struct2cell(sessionParams), 'UniformOutput', false);
else
    expParamsDef = [];
    sessionParamsDef = [];
end
noSessionParams = 1;
while noSessionParams
    dlgtitle = 'Enter session parameters';
    prompt = {'Fixation X:',...
        'Fixation Y:',...
        'Field size:',...
        'PPD X:',...
        'PPD Y:'};
    
    fieldsize(:,1) = ones(size(prompt));
    fieldsize(:,2) = 45*ones(size(prompt));
    if ~isempty(sessionParamsDef)
        answer = inputdlg(prompt, dlgtitle, fieldsize, sessionParamsDef(2:end));
    else
        answer = inputdlg(prompt, dlgtitle, fieldsize);
    end
    
    if ~isempty(answer)
        noSessionParams = 0;
        
    else
        quitAnswer = questdlg('Quit?',...
            'Stop and exit',...
            'Yes', 'No', 'No');
        if strcmp(quitAnswer, 'No')
        else
            error('User quit');
        end
    end
end

sessionParams.initials = subjectId;
sessionParams.fixationX = str2num(answer{strcmp(prompt, 'Fixation X:')});
sessionParams.fixationY = str2num(answer{strcmp(prompt, 'Fixation Y:')});
sessionParams.fieldSize = str2num(answer{strcmp(prompt, 'Field size:')});
sessionParams.PPD_X = str2num(answer{strcmp(prompt, 'PPD X:')});
sessionParams.PPD_Y = str2num(answer{strcmp(prompt, 'PPD Y:')});

%Ask user whether they want to load previous  experimental parameters

noExpParams = 1;

while noExpParams
    loadPrev = questdlg('Load previous parameters?',...
        'Experimental Parameters',...
        'Yes', 'No', 'No');
    
    switch loadPrev
        case 'No'
            
            while noExpParams == 1
                
                dlgtitle = 'Enter experimental parameters';
                prompt = {'Gain:',...
                    'Stimulus diameter (pixels):',...
                    'Stimulus duration (frames):',...
                    'No. of staircases:',...
                    'Trials per staircase:',...
                    'Red TCA X:',...
                    'Red TCA Y:',...
                    'Green TCA X:',...
                    'Green TCA Y:',...
                    'Target X:',...
                    'Target Y:',...
                    'Grid width:',...
                    'Grid height:'};
                
                fieldsize2(:,1) = ones(size(prompt));
                fieldsize2(:,2) = 45*ones(size(prompt));
                
                if ~isempty(sessionParamsDef)
                    answer2 = inputdlg(prompt, dlgtitle, fieldsize2,expParamsDef);
                else
                    answer2 = inputdlg(prompt, dlgtitle, fieldsize2);
                end
                
                if ~isempty(answer2)
                    noExpParams = 0;
                else
                    quitAnswer = questdlg('Quit?',...
                        'Stop and exit',...
                        'Yes', 'No', 'No');
                    if strcmp(quitAnswer, 'No')
                    else
                        error('User quit');
                    end
                end
            end %while
            expParams.gain = str2num(answer2{strcmp(prompt, 'Gain:')});
            expParams.stimDiam = str2num(answer2{strcmp(prompt, 'Stimulus diameter (pixels):')});
            expParams.stimDur = str2num(answer2{strcmp(prompt, 'Stimulus duration (frames):')});
            expParams.numStair = str2num(answer2{strcmp(prompt, 'No. of staircases:')});
            expParams.trialsPerStair = str2num(answer2{strcmp(prompt, 'Trials per staircase:')});
            expParams.redTCA_X = str2num(answer2{strcmp(prompt, 'Red TCA X:')});
            expParams.redTCA_Y = str2num(answer2{strcmp(prompt, 'Red TCA Y:')});
            expParams.greenTCA_X = str2num(answer2{strcmp(prompt, 'Green TCA X:')});
            expParams.greenTCA_Y = str2num(answer2{strcmp(prompt, 'Green TCA Y:')});
            expParams.targetX = str2num(answer2{strcmp(prompt, 'Target X:')});
            expParams.targetY = str2num(answer2{strcmp(prompt, 'Target Y:')});
            expParams.gridWidth = str2num(answer2{strcmp(prompt, 'Grid width:')});
            expParams.gridHeight = str2num(answer2{strcmp(prompt, 'Grid height:')});
            
        case 'Yes'
            [prevParamFile,prevParamPath,~] = uigetfile([subjectPath(1:end-1) '*.mat'], 'Choose file');
            if isstr(prevParamFile) && isstr(prevParamPath)
                try
                    load(fullfile(prevParamPath,prevParamFile), 'expParams');
                    noExpParams = 0;
                catch
                    error('Invalid file selected');
                end
            else
            end
    end
end

fileName = [sessionParams.initials '_RedGreenThresh_' datestr(now, 'mm_dd_yy_HH_MM_SS')];
save([subjectPath filesep fileName '.mat'], 'sessionParams', 'expParams')

%%
function createStim(stimsize, stimshape, intensities, extension)
% Create stimuli. Default will create a zero stimulus for frame2, cross
% for frame3 and 100% increment for frame4.
%
% USAGE
% createStimulus(stimsize, stimshape, powers)
%
% INPUT
% stimsize:     in pixels
% stimshape:    char. currently supports square or circle
% intensity:    stimulus intensity. values should be between 0
%               and 1. powers can be a single value or an array of
%               values, in which case a bmp file will be created for
%               each power.
% extension:    bmp or buf. default is bmp. bit depth is higher with
%               buf.
%
% OUTPUT
% saves bmp files into tempStimulus directory
%
if nargin < 3 || isempty(intensities)
    intensities = 1;
end
if nargin < 4
    extension = 'bmp';
end

stimdir = fullfile(pwd, 'tempStimulus');

% cycle through powers and create stimuli for each
frameN = 4;
for p = 1:length(intensities)
    intensity = intensities(p);
    if strcmp(stimshape, 'square')
        stim_im = zeros(stimsize, stimsize);
        stim_im(1:end,1:end) = intensity;
        stim_im = padarray(stim_im, [1 1],0, 'both');
    elseif strcmp(stimshape, 'circle')
        xp = -fix(stimsize / 2):fix(stimsize / 2);
        [x, y] = meshgrid(xp);
        stim_im = (x .^ 2 + y .^ 2) <= (round(stimsize / 2)) .^ 2;
        stim_im = stim_im .* intensity;
        
    end
    % write to file
    if strcmpi(extension, 'bmp')
        imwrite(stim_im, fullfile(stimdir, ...
            ['frame' num2str(frameN) '.bmp']));
    elseif strcmpi(extension, 'buf')
        stim.write_to_buf_file(stim_im, num2str(frameN), stimdir, ...
            'frame')
    end
    frameN = frameN + 1;
end

if isdir(fullfile(pwd, 'tempStimulus')) == 0
    mkdir(fullfile(pwd, 'tempStimulus'));
end

% make a blank (zero value) stimulus that is at least as large as the
% stimulus.
blank_im = stim_im > 0;

%Make it an increment
%blank_im = blank_im == 1;

% write to file
if strcmpi(extension, 'bmp')
    imwrite(blank_im, fullfile(stimdir, 'frame2.bmp'));
elseif strcmpi(extension, 'buf')
    stim.write_to_buf_file(blank_im, '2', stimdir, ...
        'frame')
end

% Make cross to record stimulus location
cross_im = stim.create_cross_img(21, 5, true);
% write to file
if strcmpi(extension, 'bmp')
    imwrite(cross_im, fullfile(stimdir, 'frame3.bmp'));
elseif strcmpi(extension, 'buf')
    stim.write_to_buf_file(cross_im, '3', stimdir, ...
        'frame')
end

% Make increment with intensity 1

increment = blank_im == 1;

%Pad again with zeros and make stimulus bigger
increment = padarray(increment, [1 1], 'both');
[i,j] = find(increment);
increment(min(i)-2:max(i)+2, min(j)-2:max(j)+2) = 1;

if strcmpi(extension, 'bmp')
    fName = ['frame' num2str(frameN) '.bmp'];
    imwrite(increment, fullfile(stimdir, fName));
elseif strcmpi(extension, 'buf')
    stim.write_to_buf_file(increment, num2str(frameN), stimdir, ...
        'frame')
end