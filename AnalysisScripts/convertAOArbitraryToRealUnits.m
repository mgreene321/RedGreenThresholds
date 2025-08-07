function valueInRealUnits = convertAOArbitraryToRealUnits(valueInArbitraryUnits,...
    maxLaserPowerWatts, lambdaMax, relativeSpectrum, relativeSpectrumS, scanAngleDeg)
% [S, valueInRealUnits] = convertAOArbitraryToRealUnits(valueInArbitraryUnits, maxLaserPowerWatts, lambdaMax, FWHM, scanAngleDeg, photoreceptors)
% 
% Convert AOSLO laser intensity in arbitrary units to real radiometric
% units as well as to isomerizations/receptor/s
%
% Input:
%    valueInArbitraryUnits: e.g., measured threshold in terms of the 0-1
%    AOM power output
%
%    maxLaserPowerWatts: laser power at cornea (measured at peak wavelength at maximum intensity)in watts
%
%    FWHM: full width at half maximum of laser spectrum
%
%    scanAngleDeg: angular subtense of scanning raster in deg
%
%    photorecptors: structure containing the following fields:
%
%                   {'species'                 }
%                   {'types'                   }
%                   {'nomogram'                }
%                   {'OSlength'                }
%                   {'ISdiameter'              }
%                   {'specificDensity'         }
%                   {'quantalEfficiency'       }
%                   {'fieldSizeDegrees'        }
%                   {'ageInYears'              }
%                   {'pupilDiameter'           }
%                   {'eyeLengthMM'             }
%                   {'lensDensity'             }
%                   {'macularPigmentDensity'   }
%                   {'axialDensity'            }
%                   {'absorbance'              }
%                   {'absorptance'             }
%                   {'preReceptoral'           }
%                   {'effectiveAbsorptance'    }
%                   {'isomerizationAbsorptance'}
%                   {'energyFundamentals'      }
%                   {'quantalFundamentals'     }
%
%               Can be created as follows:
%                 photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%                 photoreceptors = FillInPhotoreceptors(photoreceptors);
%
% Output:
%    S: [start delta num] decription of list of wavelengths
%
%    valueInRealUnits: structure containing the following fields:
%
%                   {'wattsPerLambdaIntoEye'                  }
%                   {'rawCornIrradianceWattsPerLambdaPerCm2In'}
%                   {'radianceWattsPerLambdaPerCm2Sr'         }
%                   {'retIrradianceWattsPerLambdaPerUm2'      }
%                   {'retIrradianceQuantaPerLambdaPerUm2Sec'  }
%                   {'transmittedQuantaPerLambdaPerUm2Sec'    }
%                   {'absorbedQuantaPerLambdaPerUm2Sec'       }
%                   {'effectiveQuantaPerLambdaPerUm2Sec'      }
%                   {'isomzeriationsSec'                      }
%                   {'fractionBleachedFromIsom'               }
%
% 4/16/2024 mjg Wrote it. 

%% Constants
h = 6.62607015e-34; % Planck's constant
c = 2.99792458e8; % Speed of light
pupilDiamMm = 6.5;
eyeLengthMm = 16.7; %17;
ISdiameter = 5.5; %microns, for ~2.5 deg eccentricity (e.g. Scoles et al 2014)

%% Utility calculations
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
pupilAreaCm2 = pupilAreaMm2*(10^-2);
pupilAreaM2 = pupilAreaMm2*(10^-6);
eyeLengthCm = eyeLengthMm*(10^-1);
eyeLengthM = eyeLengthMm*(10^-3);
scanAngleDeg2 = scanAngleDeg^2;
scanSolidAngle = 2*pi*(1-cosd(scanAngleDeg/2));

%% Compute power spectrum of laser light
% Define wavelength support for radiometric quantities
%S = [380 1 401];
wvl = transpose(SToWls(relativeSpectrumS));

%maxPowerPerLambda = maxLaserPowerWatts .* maxPowerPerLambda./max(maxPowerPerLambda(:));
%maxPowerPerLambda = maxLaserPowerWatts .* maxPowerPerLambda./sum(maxPowerPerLambda(:)); %sum(maxPowerPerLambda(:)) = maxPowerWatts;

responsivityFile = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\power_meter_response_curve.xlsx';

relativeSpectrum = relativeSpectrum./max(relativeSpectrum(:)); % normalize to 1
%relativeSpectrum = transpose(relativeSpectrum);
if isfile(responsivityFile)
    respTbl = readtable(responsivityFile);
    respWvl = respTbl{:,1};
    respVal = respTbl{:,2};
    responsivity = interp1(respWvl, respVal, wvl, 'linear');  
    responsivity_fun = @(lambda) interp1(respWvl, respVal, lambda, 'linear', 'extrap');
   % spectrumScalar = (maxLaserPowerWatts .* responsivity(wvl == lambdaMax))./(sum(relativeSpectrum.*responsivity, 'all', 'omitnan'));
    spectrumScalar = (maxLaserPowerWatts .* responsivity_fun(lambdaMax))./(sum(relativeSpectrum.*responsivity, 'all', 'omitnan'));
else
    spectrumScalar = maxLaserPowerWatts./sum(relativeSpectrum(:));
    warning('No responsivity data found')    
end

maxPowerPerLambda = spectrumScalar .* relativeSpectrum;
%maxPowerPerLambda = maxLaserPowerWatts .* relativeSpectrum';

% scale spectrum by threshold (watts at lambdaMax)
wattsPerLambdaIntoEye = valueInArbitraryUnits .* maxPowerPerLambda;

%% Radiometric calculations

% Following Domdei et al 2018

retinalAreaM2 = (2*tand(scanAngleDeg/2)*eyeLengthM)^2;
retIrradianceWattsPerM2lambda = wattsPerLambdaIntoEye./retinalAreaM2;
radianceWattsPerM2Srlambda = retIrradianceWattsPerM2lambda * (eyeLengthM^2)/pupilAreaM2;

% Brainard way (should yield similar results)

%  Corneal irradiance [W m^-2 lambda^-1] = power [W lambda^-1] / pupil area [m^2]
cornIrradianceWattsPerM2lambda = wattsPerLambdaIntoEye./pupilAreaM2;
% Radiance [W m^-2 Sr^-1 lambda^-1] = corneal irradiance [W m^-2 lambda^-1] / scan angle^2 [Sr]
radianceWattsPerM2Srlambda_B = cornIrradianceWattsPerM2lambda/(deg2rad(sqrt(scanAngleDeg2)))^2;
%Retinal irradiance
retIrradianceWattsPerM2lambda_B = radianceWattsPerM2Srlambda_B * pupilAreaM2/(eyeLengthM^2); 

% Check to make sure both methods agree reasonably well
relativeErrorPct = 100 * mean((retIrradianceWattsPerM2lambda_B - retIrradianceWattsPerM2lambda)/retIrradianceWattsPerM2lambda, 'all', 'omitnan');
if relativeErrorPct > 1
    warning('Relative error between Domdei and Brainard calculations is %d percent', relativeErrorPct) 
end

% Energy units
valueInRealUnits.wattsPerLambdaIntoEye = wattsPerLambdaIntoEye;
valueInRealUnits.cornIrradianceWattsPerM2lambda = cornIrradianceWattsPerM2lambda;
valueInRealUnits.radianceWattsPerM2Srlambda = radianceWattsPerM2Srlambda;
valueInRealUnits.retIrradianceWattsPerM2lambda = retIrradianceWattsPerM2lambda;

% Quantal units
energyToQuanta = (wvl * 1e-9)./(h*c);

valueInRealUnits.quantaPerSecLambdaIntoEye = wattsPerLambdaIntoEye .* energyToQuanta;
valueInRealUnits.cornIrradianceQuantaPerSecM2lambda = cornIrradianceWattsPerM2lambda .* energyToQuanta;
valueInRealUnits.radianceQautnaPerSecM2Srlambda = radianceWattsPerM2Srlambda .* energyToQuanta;
valueInRealUnits.retIrradianceQuantaPerSecM2lambda = retIrradianceWattsPerM2lambda .* energyToQuanta;

end