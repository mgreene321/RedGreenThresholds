function [paramsValues,thresh,LL,exitFlag,bootstrap]=fitPSF_RG(stimLevels,numSeen,outOfNum,PF,thresholdPercent,paramsFree,guess,lapse,niter,parametric)

% This function was written by AE Boehm to perform psychometric function fitting with the Palamedes toolbox.
% Make sure toolbox is installed and in the matlab path before using this
% function. (http://www.palamedestoolbox.org/download.html)
%
% Input parameters
%   'stimLevels': vector containing stimulus levels used.
%
%   'numSeen': vector containing for each of the entries of 'stimLevels' the 
%       number of trials a positive response (e.g., 'yes' or 'correct') was
%       given.
%
%   'OutOfNum': vector containing for each of the entries of 'stimLevels' 
%       the total number of trials.
%
%   'PF': psychometric function to be fitted. Passed as an inline function.
%       Options include:    
%           @PAL_Logistic
%           @PAL_Weibull
%           @PAL_Gumbel (i.e., log-Weibull)
%           @PAL_Quick
%           @PAL_logQuick
%           @PAL_CumulativeNormal (DEFAULT)
%           @PAL_Gumbel
%           @PAL_HyperbolicSecant
%
%   thresholdPercent: percent correct for threshold determination
%   (default = 50)
%
%   'paramsFree': 1x4 vector coding which of the four parameters of the PF 
%       [threshold slope guess-rate lapse-rate] are free parameters and 
%       which are fixed parameters (1: free, 0: fixed). 
%       (default = [0 0 1 1])
% 
%  'guess': guess rate to be used when fixed parameter. (default = 0)
%
%  'lapse': lapse rate to be used when fixed parameter. (default = 0)
%
%   'niter': number of bootstrap simulations to perform. Must be greater
%       than 0 for goodness of fit and standard errors (default = 0)
%
%   'parametric': 1 (parametric bootstrap) or 0 (nonparametric bootstrap)
%       (default = 1)
%
% Output
%   'paramsValues' The parameter values from psychometric function fit 
%       [threshold slope guess-rate lapse-rate], help PAL_PFML_Fit for details)
%
%   'thresh', Intensity at the specific thresholdPercent
%
%   'goodnessOfFit', structure containing outputs from goodness of fit
%       help PAL_PFML_GoodnessOfFit for details
%
%   'bootstrap', structure containing output from bootstrap 
%       help PAL_PFML_BootstrapParametric or
%       PAL_PFML_BootstrapNonParametric for details

% set defaults for optional params if unspecified.
if ~exist('PF','var')
    PF = @PAL_CumulativeNormal;
end
if ~exist('thresholdPercent','var')
    thresholdPercent = 50;
end
if ~exist('paramsFree','var')
    paramsFree = [1 1 0 0];
end
if ~exist('guess','var')
    guess = 0;
end
if ~exist('lapse','var')
    lapse = 0;
end
if ~exist('niter','var')
    niter = 0;
end
if ~exist('parametric','var')
    parametric = 1;
end

if(size(stimLevels,1) > size(stimLevels,2))
    stimLevels = transpose(stimLevels);
end
if(size(numSeen,1) > size(numSeen,2))
    numSeen = transpose(numSeen);
end
if(size(outOfNum,1) > size(outOfNum,2))
    outOfNum = transpose(outOfNum);
end

searchGrid.alpha = 0.01:.01:1;
searchGrid.beta = logspace(-10^-10,3,101);


lapseLimits=[0 0.5];
guessLimits=lapseLimits;

searchGrid.gamma = guess;
searchGrid.lambda = lapse;%0.05:0.01:0.5;%0.02;
    
gammaEQlambda=0;
% Do the fit


[paramsValues, LL, exitFlag] = PAL_PFML_Fit(stimLevels, numSeen, outOfNum, ...
    searchGrid, paramsFree, PF,'lapseLimits',lapseLimits,'guessLimits',guessLimits,'lapseFits', 'jAPLE');
% [paramsValues, LL, exitFlag] = PAL_PFML_Fit(stimLevels, numSeen, outOfNum, ...
%     searchGrid, paramsFree, PF,'gammaEQlambda',gammaEQlambda,'lapseLimits',[0 0.25]);

% Extract threshold at specified percent seen
thresh = PF(paramsValues, thresholdPercent./100, 'inverse');

% Do a bootstrap 
if niter > 0
    [goodnessOfFit.Dev goodnessOfFit.pDev] = PAL_PFML_GoodnessOfFit(stimLevels, numSeen, outOfNum, ...
    paramsValues, paramsFree, niter, PF, 'searchGrid', searchGrid,'lapseLimits',lapseLimits,'guessLimits',guessLimits);
    if parametric == 1
        [bootstrap.SD bootstrap.paramsSim bootstrap.LLSim bootstrap.converged] = ...
        PAL_PFML_BootstrapParametric(stimLevels, numSeen, ...
        paramsValues, paramsFree, niter, PF, ...
        'searchGrid', searchGrid,'lapseLimits',lapseLimits,'guessLimits',guessLimits);
    else
        [bootstrap.SD bootstrap.paramsSim bootstrap.LLSim bootstrap.converged] = PAL_PFML_BootstrapNonParametric(...
       stimLevels, numSeen, outOfNum, [], paramsFree, niter, PF,...
        'searchGrid',searchGrid,'lapseLimits',lapseLimits,'guessLimits',guessLimits);
    end
else
    goodnessOfFit = [];
    bootstrap.SD = [0 0];
end
end