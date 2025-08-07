function photoreceptors = photoreceptorSensitivitiesFromScratch(params)
%% Constants
loadConstants;

%%

% start with nomogram
S = [380 1 401];
S_absorbance = StockmanSharpeNomogram(S, params.nomogram.lambdaMax(3));
M_absorbance = StockmanSharpeNomogram(S, params.nomogram.lambdaMax(2));
L_absorbance = StockmanSharpeNomogram(S, params.nomogram.lambdaMax(1));
absorbance = [L_absorbance; M_absorbance; S_absorbance];

% pre-retinal filtering
lensTransmit = LensTransmittance(S, 'Human', 'CIE', params.ageInYears, params.pupilDiameter.value);
mpod = computePeakMacularPigmentDensity(params.eccentricity, params.stimSize); % for exp: 2, 2.25/60
load den_mac_bone;
macDensity = SplineSrf(S_mac_bone,den_mac_bone,S,2)';
macDensity = mpod*macDensity./max(macDensity);
macDensity(macDensity < 0) = 0;
macTransmit = 10.^(-macDensity);
preReceptoralTransmit = lensTransmit.*macTransmit;

% axialDensity
specificDensity = PhotopigmentSpecificDensity('FovealLCone', 'Human', 'Rodieck');
axialDensity = specificDensity.*params.OSlength.value.*ones(3,1);

% compute absorptance from absorbance
absorptance = 1 - 10.^(-diag(axialDensity)*absorbance);
effectiveAbsorptance = absorptance .* (ones(size(absorptance,1),1)*preReceptoralTransmit);

quantalEfficiency = [0.667 0.667 0.667]';

for i = 1:size(effectiveAbsorptance,1)
    isomerizationAbsorptance(i,:) = quantalEfficiency(i) * ...
        effectiveAbsorptance(i,:);
end


energyFundamentals = EnergyToQuanta(S,isomerizationAbsorptance')';
mx = max(energyFundamentals,[],2);
energyFundamentals = diag(1./mx)*energyFundamentals;

%% Compute normalized quantal sensitivities (aka cone fundamentals in quantal units)
quantalFundamentals = isomerizationAbsorptance;
mx = max(quantalFundamentals,[],2);
quantalFundamentals = diag(1./mx)*quantalFundamentals;

%% Account for bleaching by projector background

coneAperture = coneAperture./max(coneAperture(:));
projectorRGB = whiteRGB;
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix= pixelSideLengthUm.^2;

rgb_spectra = cal.rgb_spectra;
wavelength_sampling = transpose(cal.wavelength_sampling);
norm_lum_output = cal.norm_lum_output;

R_max = rgb_spectra(norm_lum_output(:,1) == 1,:,1);
G_max = rgb_spectra(norm_lum_output(:,2) == 1,:,2);
B_max = rgb_spectra(norm_lum_output(:,3) == 1,:,3);

R_max = R_max(1,:);
G_max = G_max(1,:);
B_max = B_max(1,:);

projectorRGB_norm = projectorRGB./255;

R_scalar = projectorRGB_norm(1);
G_scalar = projectorRGB_norm(2);
B_scalar = projectorRGB_norm(3);
% projector spd
projectorSpectrum = transpose(R_scalar.*R_max + G_scalar.*G_max  + B_scalar.*B_max);


%projectorSpectrum = SplineSpd(wavelength_sampling', projectorSpectrum', wvl);

% spline fundamentals to wavelength domain

isomerizationAbsorptanceSplined = SplineCmf(S, isomerizationAbsorptance, wavelength_sampling,2);

%[projectorRadianceWattsPerM2Srlambda, ~] = LumToRadiance(projectorSpectrum, S, projectorLum); %now, 683*dot(T_vLambda, projectorRadianceWattsPErM2Srlambda) = projectorLum
[projectorRadianceWattsPerM2Srlambda, ~] = LumToRadiance(projectorSpectrum, WlsToS(wavelength_sampling), projectorLum); %now, 683*dot(T_vLambda, projectorRadianceWattsPErM2Srlambda) = projectorLum

% retinal irradiance = radiance * pupilArea/eyeLength^2
%projectorRetIrradianceWattsPerM2lambda = RadianceAndPupilAreaEyeLengthToRetIrradiance(projectorRadianceWattsPerM2Srlambda, S, pupilAreaCm2,eyeLengthCm);
projectorRetIrradianceWattsPerM2lambda = RadianceAndPupilAreaEyeLengthToRetIrradiance(projectorRadianceWattsPerM2Srlambda, WlsToS(wavelength_sampling), pupilAreaCm2,eyeLengthCm);

projectorRetIrradianceQuantaPerSecM2lambda = projectorRetIrradianceWattsPerM2lambda .* (wavelength_sampling .* 1e-9)/(h*c);
projectorRetIrradianceQuantaPerSecUm2lambda = projectorRetIrradianceQuantaPerSecM2lambda *1e-12;

% quanta of each wavelength arriving per second within AO pixel sized area on retina
projectorRetIrradianceQuantaPerSecPixlambda = projectorRetIrradianceQuantaPerSecUm2lambda.* Um2PerPix;

for i = 1:numel(wavelength_sampling)
    transmittedQuantaPerSecLambda(i) = sum(coneAperture.*projectorRetIrradianceQuantaPerSecPixlambda(i), 'all');
    %transmittedQuantaPerSecLambda(i) = sum(coneAreaUm2*projectorRetIrradianceQuantaPerSecUm2lambda(i), 'all');
end

for rr = 1:3
    isomerizationsSec(rr,:) = dot(transmittedQuantaPerSecLambda, isomerizationAbsorptanceSplined(rr,:));
    fractionBleachedFromIsom(rr,:) = ComputePhotopigmentBleaching(isomerizationsSec(rr,:),'cones','isomerizations','Boynton');
end

axialDensityBleached = axialDensity .* (1-fractionBleachedFromIsom);

absorptanceBleached = 1 - 10.^(-diag(axialDensityBleached)*absorbance);
effectiveAbsorptanceBleached = absorptanceBleached .* (ones(size(absorptanceBleached,1),1)*preReceptoralTransmit);

quantalEfficiency = [0.667 0.667 0.667]';

for i = 1:size(effectiveAbsorptance,1)
    isomerizationAbsorptanceBleached(i,:) = quantalEfficiency(i) * ...
        effectiveAbsorptanceBleached(i,:);
end

energyFundamentalsBleached = EnergyToQuanta(S,isomerizationAbsorptanceBleached')';
mx = max(energyFundamentalsBleached,[],2);
energyFundamentalsBleached = diag(1./mx)*energyFundamentalsBleached;

%% Compute normalized quantal sensitivities (aka cone fundamentals in quantal units)
quantalFundamentalsBleached = isomerizationAbsorptanceBleached;
mx = max(quantalFundamentalsBleached,[],2);
quantalFundamentalsBleached = diag(1./mx)*quantalFundamentalsBleached;

photoreceptors.isomerizationAbsorptance = isomerizationAbsorptance;
photoreceptors.isomerizationAbsorptanceBleached = isomerizationAbsorptanceBleached;
photoreceptors.projectorRetIrradianceQuantaPerSecUm2lambda = projectorRetIrradianceQuantaPerSecUm2lambda;
photoreceptors.nomogram.S = S;

end