clear; close all;


% Macular pigment peak optical density as a function of field size (CIE
% 16x:2005)
Dt = @(fs) 0.485*exp(-fs./6.132);

% Sample field sizes (deg) and corresponding radii
fieldSizes = 0:0.5:10;
radii = fieldSizes./2;

% Compute optical density in the region between two fields of radii R1 and
% R2, where R2 > R1

for i = 1:numel(radii)-1

    % radii of fields
    R1 = radii(i);
    R2 = radii(i+1);

    Ravg(i) = mean([R1 R2]);

    % peak optical densities
    Dt1 = Dt(fieldSizes(i)); 
    Dt2 = Dt(fieldSizes(i+1));

    % transmittance within region between fields
    T(i) = ((R2^2)*10^-Dt2 - (R1^2)*10^-Dt1)./(R2^2 - R1^2);

    % optical density = -log10(transmittance)
    OD(i) = -log10(T(i));

end

% plot OD vs. eccentricity
midPoints = radii(1:end-1) + 0.25/2;
figure, hold on
plot(Ravg, OD, 'ko', 'LineStyle', 'none', 'LineWidth', 2)
xlabel('Retinal eccentricity (deg)'); ylabel('Macular pigment optical density')
set(gca, 'fontSize', 14)

%Hammond_MPOD = importdata('C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Hammond_MPOD.mat');

% compare with data from Hammond et al. (1997, JOSA A)
Hammond_MPOD = [0.0011    0.5811
                0.5222    0.3913
                1.0043    0.2913
                1.9908    0.1260
                3.0026    0.0796
                4.0137    0.0143
                5.5017   -0.0010];

plot(Hammond_MPOD(:,1), Hammond_MPOD(:,2), 'LineWidth', 2)
legend {CIE formula} {Hammond et al. 1997}