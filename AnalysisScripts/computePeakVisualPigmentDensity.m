function vpod = computePeakVisualPigmentDensity(eccentricity, annulusWidth)
Dt_LM = @(fs) 0.38 + 0.54*exp(-fs/1.333);
Dt_S = @(fs) 0.3+0.45*exp(-fs/1.333);

R1 = eccentricity - annulusWidth;
R2 = eccentricity + annulusWidth;

Dt_LM1 = Dt_LM(2*R1);
Dt_LM2 = Dt_LM(2*R2);

Dt_S1 = Dt_S(2*R1);
Dt_S2 = Dt_S(2*R2);

T_LM = ((R2^2)*10^-Dt_LM2 - (R1^2)*10^-Dt_LM1)./(R2^2 - R1^2);
vpod_LM = -log10(T_LM);

T_S = ((R2^2)*10^-Dt_S2 - (R1^2)*10^-Dt_S1)./(R2^2 - R1^2);
vpod_S = -log10(T_S);

vpod(1:2) = vpod_LM;
vpod(3) = vpod_S;

vpod = vpod';

end