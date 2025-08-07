function params = getSubjectParams(subjectId)

params.macularPigmentDensity.source = 'CIE';
params.lensDensity.source = 'CIE';
params.stimSize = 2.25/60;
params.fractionPigmentBleached.value = zeros(3,1);
params.nomogram.lambdaMax = [558.9 530.3 420.7]';
params.pupilDiameter.value = 6.5; % mm


switch subjectId
    case '10001R'
        params.ageInYears = 40.5;
        params.eccentricity = 2;
        params.OSlength.value = 35.7;
        params.propL = 0.5954;
    case '20217R'
        params.ageInYears = 27;
        params.eccentricity = 2.5;
        params.OSlength.value = 32.9;
        params.propL = 0.7852;
end

end