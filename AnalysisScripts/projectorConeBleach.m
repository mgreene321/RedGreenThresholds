function [fractionBleachedFromIsom, isomerizationsSec, projectorRetIrradianceQuantaPerSecUm2lambda] = projectorConeBleach(projectorRGB, projectorLum, projectorParams, eyeParams, photoreceptors)

% projectorParams: RGB, luminance, calibration, ND
% eyeParams: pupilAreaCm2, eyeLengthCm, ISareaUm2

h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light

S = photoreceptors.nomogram.S;
wvl = SToWls(S);

rgb_spectra = projectorParams.cal.rgb_spectra;
wavelength_sampling = transpose(projectorParams.cal.wavelength_sampling);
norm_lum_output = projectorParams.cal.norm_lum_output;

%pupilAreaCm2 = eyeParams.pupilAreaCm2;
pupilDiamMm = 6.5;
pupilDiamCm = pupilDiamMm*(10^-1);
pupilAreaCm2 = pi *(pupilDiamCm/2)^2;
eyeLengthCm = eyeParams.eyeLengthCm;
eyeLengthMm = eyeLengthCm*10;
%ISdiameterUm = eyeParams.ISdiameterUm;
coneAperture = eyeParams.coneAperture;
coneAperture =  coneAperture./max(coneAperture(:));

% incorrect, before 6/30/25
% projectorSpectrum = rgb_spectra(projectorRGB(1)+1,:,1) + ...
%                     rgb_spectra(projectorRGB(2)+1,:,2) + ...
%                     rgb_spectra(projectorRGB(3)+1,:,3);
                
%spline spectrum to wavelength representation of photoreceptors

R_max = rgb_spectra(norm_lum_output(:,1) == 1,:,1);
G_max = rgb_spectra(norm_lum_output(:,2) == 1,:,2);
B_max = rgb_spectra(norm_lum_output(:,3) == 1,:,3);
% 
% R_max = cal.rgb_spectra(end,:,1);
% G_max = cal.rgb_spectra(end,:,2);
% B_max = cal.rgb_spectra(end,:,3);

% get rid of duplicates
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

photoreceptors.isomerizationAbsorptance = SplineCmf(S, photoreceptors.isomerizationAbsorptance, wavelength_sampling,2);

%[projectorRadianceWattsPerM2Srlambda, ~] = LumToRadiance(projectorSpectrum, S, projectorLum); %now, 683*dot(T_vLambda, projectorRadianceWattsPErM2Srlambda) = projectorLum
[projectorRadianceWattsPerM2Srlambda, ~] = LumToRadiance(projectorSpectrum, WlsToS(wavelength_sampling), projectorLum); %now, 683*dot(T_vLambda, projectorRadianceWattsPErM2Srlambda) = projectorLum


% retinal irradiance = radiance * pupilArea/eyeLength^2
%projectorRetIrradianceWattsPerM2lambda = RadianceAndPupilAreaEyeLengthToRetIrradiance(projectorRadianceWattsPerM2Srlambda, S, pupilAreaCm2,eyeLengthCm);
projectorRetIrradianceWattsPerM2lambda = RadianceAndPupilAreaEyeLengthToRetIrradiance(projectorRadianceWattsPerM2Srlambda, WlsToS(wavelength_sampling), pupilAreaCm2,eyeLengthCm);

projectorRetIrradianceQuantaPerSecM2lambda = projectorRetIrradianceWattsPerM2lambda .* (wavelength_sampling .* 1e-9)/(h*c);
projectorRetIrradianceQuantaPerSecUm2lambda = projectorRetIrradianceQuantaPerSecM2lambda *1e-12; 

% convert to pixel

ppd = 560;
pixelSideLengthDeg = 1/ppd;
pixelSideLengthMm = DegreesToRetinalMM(pixelSideLengthDeg, eyeLengthMm);
pixelSideLengthUm = pixelSideLengthMm*1e3;
Um2PerPix = pixelSideLengthUm.^2;

% quanta of each wavelength arriving per second within AO pixel sized area on retina
projectorRetIrradianceQuantaPerSecPixlambda = projectorRetIrradianceQuantaPerSecUm2lambda.* Um2PerPix; 
% scalar = sum(projectorRetIrradianceQuantaPerSecPix(:)); %sum across wavelengths
% normProj = projectorRetIrradianceQuantaPerSecPix./scalar;

% for each wavelength, compute the quanta/s collected by the Gaussian cone
% aperture

for i = 1:numel(wavelength_sampling)
    transmittedQuantaPerSecLambda(i) = sum(coneAperture.*projectorRetIrradianceQuantaPerSecPixlambda(i), 'all');

end



L = photoreceptors.isomerizationAbsorptance(1,:);
M = photoreceptors.isomerizationAbsorptance(2,:);
S = photoreceptors.isomerizationAbsorptance(3,:);

% 
% Lproj = dot(L', normProj);
% Mproj = dot(M', normProj);
% Sproj = dot(S', normProj);

% coneProj = [Lproj Mproj Sproj];


% OLD WAY (turns out to be the same as the new way it seems)

% for rr = 1:3
%     %T_quantalIsomerizations = photoreceptors.isomerizationAbsorptance(rr,:);
%     %ISareaUm2 = pi*(ISdiameterUm/2)^2;
%     % isomerizations [1/sec/receptor] = inner segment area [um^2/receptor] * dot((quantal efficiency [1/quanta/lambda] * absorptance []), retinal irradiance [quanta/lambda/um^2/sec])
%     %isomerizationsSec(rr,:) = ISareaUm2*dot(T_quantalIsomerizations,projectorRetIrradianceQuantaPerSecUm2lambda);
%     isomerizationsSec(rr,:) = sum(coneAperture.*coneProj(rr).*scalar, 'all');
%     fractionBleachedFromIsom(rr,:) = ComputePhotopigmentBleaching(isomerizationsSec(rr,:),'cones','isomerizations','Boynton');
% end

for rr = 1:3
    isomerizationsSec(rr,:) = dot(transmittedQuantaPerSecLambda, photoreceptors.isomerizationAbsorptance(rr,:));
    fractionBleachedFromIsom(rr,:) = ComputePhotopigmentBleaching(isomerizationsSec(rr,:),'cones','isomerizations','Boynton');
end
% 
% photoreceptorsBleached = DefaultPhotoreceptors('LivingHumanFovea');
% photoreceptorsBleached.nomogram.lambdaMax = photoreceptors.nomogram.lambdaMax;
% 
% photoreceptorsBleached.fractionPigmentBleached.value = fractionBleachedFromIsom;
% 
% photoreceptorsBleached = FillInPhotoreceptors(photoreceptorsBleached);
end
