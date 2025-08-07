function scalar = computeProjectorScalar(RGB, power, lambda)
S = [380 1 401];
wvl = SToWls(S);

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
responsivityFile = 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\RedGreenThresholds\power_meter_response_curve.xlsx';

 respTbl = readtable(responsivityFile);
    respWvl = respTbl{:,1};
    respVal = respTbl{:,2};
    responsivity = interp1(respWvl, respVal, wvl);  

relativeSpectrum = cal.rgb_spectra(RGB(1)+1,:,1) + ...
                   cal.rgb_spectra(RGB(2)+1,:,2) + ...
                   cal.rgb_spectra(RGB(3)+1,:,3);

 relativeSpectrum = SplineSpd(cal.wavelength_sampling', relativeSpectrum',wvl);
               
scalar =  power.*responsivity(wvl==lambda)./(sum(relativeSpectrum.*responsivity, 'all', 'omitnan'));
end