function optimal_wL = findOptimalConeWeights

%% Constants

loadConstants;

root = getenv('USERPROFILE');
projectDir = dir([root '\Documents\**' filesep '*RedGreenThresholds']);
RedGreenThresholdsPath = fullfile(projectDir(1).folder, projectDir(1).name);
ADT = importdata(fullfile(RedGreenThresholdsPath, 'analyzedDataTable.mat'));
ADT = ADT(ADT.StimDurFrames == 3,:);
subjects = unique(ADT.SubjectID);
propL = ADT.L./(ADT.L + ADT.M);

convertThresholdsToQuanta;
coneSpectralPeaks = [555.5 530.3 420.7; 558.9 530.3 420.7; 563.4 530.3 420.7];

for subj = 1:numel(subjects)
    subjIdx = strcmpi(ADT.SubjectID, subjects{subj});
    params = getSubjectParams(subjects{subj});
    subjPropL = propL(subjIdx);
    realData{subj} = retIrradGreenQuantaPerSecM2(subjIdx)./retIrradRedQuantaPerSecM2(subjIdx);
    n = numel(realData{subj});
    for p = 1:size(coneSpectralPeaks,1)
        params.nomogram.lambdaMax = transpose(coneSpectralPeaks(p,:));
        tempFun = @(wL) RMSE(params, subjPropL, wL, realData{subj});
        x0 = 1;
        optimal_wL{subj,p} = fmincon(tempFun, x0, [],[],[],[],zeros(size(x0)),[]);
       disp(RMSE(params, subjPropL, optimal_wL{subj,p}, realData{subj}));
    end
end



function err = RMSE(params, subjPropL, wL, realData)
[RedSensitivity, GreenSensitivity] = simplePigmentModel(params, subjPropL, wL);
RG = RedSensitivity./GreenSensitivity;
n = numel(realData);
err = sqrt((1/n).*sum((RG - realData).^2));




