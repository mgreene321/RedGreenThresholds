function customPhotoreceptors = customConeSensitivities(params)

photoreceptors = DefaultPhotoreceptors('CIE2Deg');
%photoreceptors.OSlength.source = 'None';
%photoreceptors.specificDensity.source = 'None';
photoreceptors.macularPigmentDensity.source = params.macularPigmentDensity.source;
photoreceptors.lensDensity.source = params.lensDensity.source; % CIE
photoreceptors.ageInYears = params.ageInYears;
%photoreceptors.ISdiameter.value = params.ISdiameter.value;
photoreceptors.OSlength.value = params.OSlength.value;
photoreceptors.nomogram.lambdaMax = params.nomogram.lambdaMax;
photoreceptors.fractionPigmentBleached.value = params.fractionPigmentBleached.value;
photoreceptors.pupilDiameter.value = params.pupilDiameter.value;

stimSize = params.stimSize; % in deg
eccentricity = params.eccentricity;

% compute peak macular pigment optical density based on eccentricity and
% field size


S = photoreceptors.nomogram.S;

mpod = computePeakMacularPigmentDensity(eccentricity, stimSize);
%vpod = computePeakVisualPigmentDensity(eccentricity, stimSize);
%photoreceptors.axialDensity.value = vpod;

photoreceptors = FillInPhotoreceptors(photoreceptors);

customPhotoreceptors = photoreceptors;

customPhotoreceptors.macularPigmentDensity.density = mpod .* ...
    customPhotoreceptors.macularPigmentDensity.density./max(customPhotoreceptors.macularPigmentDensity.density);

customPhotoreceptors.macularPigmentDensity.density(customPhotoreceptors.macularPigmentDensity.density < 0) = 0;

customPhotoreceptors.macularPigmentDensity.transmittance = 10.^-customPhotoreceptors.macularPigmentDensity.density;

customPhotoreceptors.preReceptoral.transmittance = customPhotoreceptors.lensDensity.transmittance .* ...
    customPhotoreceptors.macularPigmentDensity.transmittance;

customPhotoreceptors.effectiveAbsorptance = customPhotoreceptors.absorptance .* ...
    (ones(size(customPhotoreceptors.absorptance,1),1)*customPhotoreceptors.preReceptoral.transmittance);

for i = 1:size(customPhotoreceptors.effectiveAbsorptance,1)
    customPhotoreceptors.isomerizationAbsorptance(i,:) = customPhotoreceptors.quantalEfficiency.value(i) * ...
        customPhotoreceptors.effectiveAbsorptance(i,:);
end

customPhotoreceptors.energyFundamentals = EnergyToQuanta(S,customPhotoreceptors.isomerizationAbsorptance')';
mx = max(customPhotoreceptors.energyFundamentals,[],2);
customPhotoreceptors.energyFundamentals = diag(1./mx)*customPhotoreceptors.energyFundamentals;

%% Compute normalized quantal sensitivities (aka cone fundamentals in quantal units)
customPhotoreceptors.quantalFundamentals = customPhotoreceptors.isomerizationAbsorptance;
mx = max(customPhotoreceptors.quantalFundamentals,[],2);
customPhotoreceptors.quantalFundamentals = diag(1./mx)*customPhotoreceptors.quantalFundamentals;

end