% load projector calibration file

cal = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\cal_03_26_2024.mat');
cal_wls = cal.wavelength_sampling';

lut = cal.lookup_table;

projectorRGB = [82 90 128];

% projector primary spds at max intensity
R_max = cal.rgb_spectra(cal.norm_lum_output(:,1) == 1,:,1);
G_max = cal.rgb_spectra(cal.norm_lum_output(:,2) == 1,:,2);
B_max = cal.rgb_spectra(cal.norm_lum_output(:,3) == 1,:,3);
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
projectorSpectrum = R_scalar.*R_max + G_scalar.*G_max  + B_scalar.*B_max;


% projectorSpectrum = rgb_spectra(cal.int_list==projectorRGB(1),:,1) + ...
%                     rgb_spectra(cal.int_list==projectorRGB(2),:,2) + ...
%                     rgb_spectra(cal.int_list==projectorRGB(3),:,3);

%% smith pokorny
% load Smith Pokorny fundamentals

smithPokorny = readmatrix('C:\Users\TutenLab-NUC3\Downloads\sp.csv');
sp_wls = smithPokorny(:,1);
smithPokorny = 10.^(smithPokorny(:,2:end));


% spline fundamentals wavelength domain

smithPokorny = SplineCmf(sp_wls, smithPokorny', cal_wls,2);

l_sp = smithPokorny(1,:);
m_sp = smithPokorny(2,:);
s_sp = smithPokorny(3,:);
s_sp = s_sp.*0.01608/0.00801;

lproj_sp = dot(projectorSpectrum, l_sp);
mproj_sp = dot(projectorSpectrum, m_sp);
sproj_sp = dot(projectorSpectrum, s_sp);

rproj_sp = lproj_sp/(lproj_sp+mproj_sp);
bproj_sp = sproj_sp/(lproj_sp+mproj_sp);
%% stockman sharpe

% load StockmanSharpe 

StockSharpe2000Struct = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Projector_stuff\StockSharpe2000_2deg.mat');
StockSharpe2000 = StockSharpe2000Struct.StockSharpe2000;


ss_wls = StockSharpe2000(:,1);
StockSharpe2000 = 10.^StockSharpe2000(:,2:end);

stockSharpe2000 = SplineCmf(ss_wls, StockSharpe2000', cal_wls,2);

l_ss = stockSharpe2000(1,:);
m_ss = stockSharpe2000(2,:);
s_ss = stockSharpe2000(3,:);

% rescale
l_ss = 0.7 .* l_ss ./ sum(l_ss);
m_ss = 0.3 .* m_ss ./ sum(m_ss);
s_ss = 1.0 .* s_ss ./ sum(s_ss);

lproj_ss = dot(projectorSpectrum, l_ss);
mproj_ss = dot(projectorSpectrum, m_ss);
sproj_ss = dot(projectorSpectrum, s_ss);

rproj_ss = lproj_ss/(lproj_ss+mproj_ss);
bproj_ss = sproj_ss/(lproj_ss+mproj_ss);

% as a check compute RGB to LMS matrix and make sure we get expected RGB
RGB2LMS = [l_ss; m_ss; s_ss] * squeeze(cal.rgb_spectra(end,:,:));
RGB_eew = inv(RGB2LMS)*[0.7;0.3;1.0];
RGB2LMS = RGB2LMS.*max(RGB_eew);

bitsScalar = 255;
projectorLLM = 0.7;
projectorSLM = 1.0;
projectorLplusM = 0.5;
rgb = round(inv(RGB2LMS) * [projectorLLM*projectorLplusM; (1-projectorLLM)*projectorLplusM; projectorSLM*projectorLplusM].*bitsScalar);