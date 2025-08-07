%Add whatever folder contains Stockman & Sharpe cone fundamentals
addpath 'C:\Users\TutenLab-NUC3\Documents\Tuten_Lab\Unique_yellow\analysistoolbox\'
load('StockSharpe2000');
StockSharpe2000(:,2:end) = 10.^StockSharpe2000(:,2:end);
wvl = StockSharpe2000(:,1);
l =StockSharpe2000(:,2);
m = StockSharpe2000(:,3);
s = StockSharpe2000(:,4);

%Interpolate cone fundamentals 
%l = @(lambda) interp1(wvl, l, lambda, 'linear');
%m = @(lambda) interp1(wvl, m, lambda, 'linear');
%% Sensitivity vs. proportion L, LUMINANCE MECHANISM

%Sensitivity if L and M-cone inputs are combined linearly, see Eq. 2 in
%Overleaf document 
S = @(lambda, propL,N) N.*(propL.*(l(wvl==lambda)) + (1-propL).*m(wvl==lambda));

%Proportion L: L/(L+M)
propL = linspace(0,1,1e3);

%Plotting------------------------------------------------------------------
figure, plot(propL, log10(S(680, propL)), 'r--', 'LineWidth', 2);
hold on; grid on;
plot(propL, log10(S(543, propL)), 'g--', 'LineWidth',2);
plot(propL, log10(S(680,propL)./S(543,propL)), 'k-', 'LineWidth',2);

%set(gca, 'YScale', 'log');
ylim([-3 0.5])
xlabel('L/(L+M)'); ylabel('Predicted log sensitivity')
legend {S_{680}} {S_{543}}
set(gca,'TickLabelInterpreter', 'latex')
set(gca, 'TickDir', 'out', 'FontSize', 12)
yLabLeft = {'-3', '2.5', '-2', '-1.5', '-1', '-0.5', '$\bar{m}(543)$', '', '0.5'};
yTicksLeft = sort([-3:0.5:0.5 log10(m(543))]);
set(gca, 'YTick', yTicksLeft, 'YTickLabel', yLabLeft);
yyaxis right
yLabRight = {'$\bar{m}(680)$', '$\bar{l}(680)$', '$\bar{l}(543)$'};
yTicksRight = sort(log10([l(680) m(680) l(543)]));
ylim([-3 0.5])
set(gca, 'YColor', [0 0 0], 'YTick', yTicksRight, 'YTickLabel', yLabRight);
set(gca, 'TickDir', 'out', 'FontSize', 12)


%% Proprtion L vs. Red/Green Sensitivity Ratio

%Minimum possible prop L is when N_L = 0, in which case S_680/S_543 =
%m(680)/m(543):

min_RGsensitivity_ratio = m(680)./m(543);

%Maximum possible propL is when N_M = 0, in which case S_680/S_543 =
%l(680)/l(543)

max_RG_sensitivity_ratio = l(680)./l(543);

RG_sensitivity_ratio = linspace(min_RGsensitivity_ratio, max_RG_sensitivity_ratio, 1e3);

S_680 = RG_sensitivity_ratio;
S_543 = ones(size(S_680));

%Solve system of equations to get number of L and M-cones, given red/green
%sensitivity ratio 
A = [l(680) m(680); l(543) m(543)];
b = [S_680; S_543];
x = A\b; %basically inv(A)*b
N_L = x(1,:); N_M = x(2,:);

propL = N_L./(N_L + N_M);

%Plotting------------------------------------------------------------------ 
figure, plot(log10(RG_sensitivity_ratio), propL, 'LineWidth', 2, 'Color', 'k')
set(gca,'TickLabelInterpreter', 'latex')
grid on;
set(gca, 'TickDir', 'out', 'FontSize', 12)
xlabel('log S_{680}/S_{543}')
ylabel('Predicted L/(L+M)');

%% Predicted psychometric functions

n = 6; %Number of quanta required to activate a cone 
quanta = linspace(0,1000, 1e4);
QL_red = zeros(size(quanta));
QM_red = QL_red;
QL_green = QL_red;
QM_green = QL_red;

for k = 0:n-1
    l_poiss_red = (exp(-l(680).*quanta).*(l(680).*quanta).^k)./factorial(k);
    m_poiss_red = (exp(-m(680).*quanta).*(m(680).*quanta).^k)./factorial(k);
    l_poiss_green = (exp(-l(543).*quanta).*(l(543).*quanta).^k)./factorial(k);
    m_poiss_green = (exp(-m(543).*quanta).*(m(543).*quanta).^k)./factorial(k);
    QL_red = QL_red + l_poiss_red;
    QM_red = QM_red + m_poiss_red;
    QL_green = QL_green + l_poiss_green;
    QM_green = QM_green + m_poiss_green;
end

%In the following analysis, set N_M = 1, so propL = N_L./(1+N_L);

PsiFunRed = @(N_L, N_M) 1 - (QL_red.^N_L).*(QM_red.^N_M);
PsiFunGreen = @(N_L, N_M) 1 - (QL_green.^N_L).*(QM_green.^N_M);

N = 100;
LM = 2.^(-1:3); %L:M ratios of interest
N_L = N.*LM; N_M = N.*ones(size(LM));

propL = N_L./(N_L+N_M);

%Plotting------------------------------------------------------------------
figure; hold on; grid on;
set(gca,'TickLabelInterpreter', 'latex')
set(gca, 'TickDir', 'out', 'FontSize', 12)

for i = 1:length(propL)
    %plot(log10(quanta), PsiFunRed(propL(i)), 'Color',[1 1-i./(length(propL)) 1-i./(length(propL))] , 'LineWidth', 1.25)
    %plot(log10(quanta), PsiFunGreen(propL(i)), 'Color',[1-i./(length(propL)) 1 1-i./(length(propL))] , 'LineWidth', 1.25)
    plot(log10(quanta), PsiFunRed(N_L(i), N_M(i)), 'Color',[1 1-i./(length(LM)) 1-i./(length(LM))] , 'LineWidth', 1.25)
    plot(log10(quanta), PsiFunGreen(N_L(i), N_M(i)), 'Color',[1-i./(length(LM)) 1 1-i./(length(LM))] , 'LineWidth', 1.25)
end
xlim([0 3])
xlabel('log relative quanta');
ylabel('P(''Seen'')')
text(1, 0.4, 'L:M ratios: 0.5, 1, 2, 4, 8', 'FontSize', 12)

%% Sensitivity vs. proportion L, CHROMATIC MECHANISM

% 
% 
% %For simplicity, assume total number of L-cones is 1, total M-cones is 1
% 
% UY = 580; %nm
% propL_global = 0.5;
% 
% k_global = l(UY)./m(UY);
% R_global = @(propL, lambda)  propL.*l(lambda) - k_global.*(1-propL).*m(lambda);
% 
% k_local = @(propL) propL.*l(UY)./((1-propL).*m(UY));
% R_local = @(propL, lambda)  propL.*l(lambda) - k_local(propL).*(1-propL).*m(lambda);
% 
% propL = linspace(0,1,1e3);
% 
% figure, plot(propL, R_global(propL, UY));
% hold on
% plot(propL, R_local(propL, UY));
% 
% S_global = @(propL, lambda) abs(R_global(propL, lambda));
% S_local = @(propL, lambda) abs(R_local(propL, lambda));
% 

