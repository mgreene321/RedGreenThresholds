function mpod = computePeakMacularPigmentDensity(eccentricity, annulusWidth)


% mpod for a foveally field size x is computePeakMacularPigmentDensity(x/4, x/4) 
Dt = @(fs) 0.485*exp(-fs./6.132);

R1 = eccentricity - annulusWidth;
R2 = eccentricity + annulusWidth;

Dt1 = Dt(2*R1);
Dt2 = Dt(2*R2);

T = ((R2^2)*10^-Dt2 - (R1^2)*10^-Dt1)./(R2^2 - R1^2);
mpod = -log10(T);

end