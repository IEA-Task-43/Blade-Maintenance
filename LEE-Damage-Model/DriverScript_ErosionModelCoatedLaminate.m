filename = 'NREL 5MW Power Curve.csv';
opts = detectImportOptions(filename);
opts.SelectedVariableNames = [1,7];
opts.DataLines = [2,24];
TurbineModel = readmatrix(filename,opts);
RotorDiameter = 126; % m
BladeLifetime = 20; % years

% filename = 'DTU 10 MW Power Curve V2.csv';
% opts = detectImportOptions(filename);
% opts.SelectedVariableNames = [1,2];
% opts.DataLines = [2,19];
% TurbineModel = readmatrix(filename,opts);
% RotorDiameter = 178.3; % m

%% Ranges of wind speed, rain intensity, and droplet size
phi_d = linspace(0.01, 6, 50);      % mm
% 
I = linspace(0.01, 50, 50);         % mm/hr

U_hub = linspace(0.0, 30, 50);      % m/s


%% Define Rainfall condition ranges - De Bilt
P_I(1).Min = 0;
P_I(1).Max = 2.5;
P_I(1).Probability = 10.29;
P_I(2).Min = 2.5;
P_I(2).Max = 10;
P_I(2).Probability = 1.35;
P_I(3).Min = 10;
P_I(3).Max = 50;
P_I(3).Probability = 0.0910;
P_I(4).Min = 50;
P_I(4).Max = NaN;
P_I(4).Probability = 0.00069;
%%

% %% Define Rainfall condition ranges - De Kooy
% P_I(1).Min = 0;
% P_I(1).Max = 2.5;
% P_I(1).Probability = 10.08;
% P_I(2).Min = 2.5;
% P_I(2).Max = 10;
% P_I(2).Probability = 1.364;
% P_I(3).Min = 10;
% P_I(3).Max = 50;
% P_I(3).Probability = 0.0801;
% P_I(4).Min = 50;
% P_I(4).Max = NaN;
% P_I(4).Probability = 0.00046;
%%


%% De Bilt
% Statistical Distributions

% Wind Distribution Variables
stat.alpha_u = 1.8763;
stat.beta_u = 5.2162;

% Rain Intensity Distribution  Variables
stat.mu_I = -0.1816;
stat.sig_I = 0.8617;

% DSD Variables
stat.N = 4.5670; 
stat.q = 0.1404;
stat.A = 0.4811; 
stat.p = 0.1186; 

% %% De Kooy
% % Statistical Distributions

% % Wind Distribution Variables
% stat.alpha_u = 1.9331;
% stat.beta_u = 8.9419;

% % Rain Intensity Distribution  Variables
% stat.mu_I = -0.1445;
% stat.sig_I = 0.8275;

% % DSD Variables
% stat.N = 2.83; 
% stat.q = -0.0953;
% stat.A = 1.03; 
% stat.p = 0.138; 

%% Material & Structural Constants

% %Substrate
% mat.rho_s     = 1900;         % density - Laminate (kg/m^3)
% mat.E_s       = 3.53e9;       % modulus of elasticity - matrix(Pa)
% mat.sig_ult_s = 76.3e6;       % Ultimate strength - matrix (Pa)
% mat.m_s       = 10.32;        % Wohler Exponent - matrix
% mat.nu_sm     = 0.347;        % Poisson's ratio - matrix
% mat.c_s       = 2740;         % speed of sound - Laminate (m/s)

%Substrate
mat.rho_s     = 1900;         % density - Laminate (kg/m^3)
mat.E_s       = 2.41e9;       % modulus of elasticity - matrix(Pa)
mat.sig_ult_s = 67.3e6;       % Ultimate strength - matrix (Pa)
mat.m_s       = 7;            % Wohler Exponent - matrix
mat.nu_sm     = 0.3; %0.399;        % Poisson's ratio - matrix
mat.c_s       = 2740;  %(mat.E22_s/mat.rho_s)^0.5;          % speed of sound - Laminate (m/s)

%Coating
mat.rho_c     = 1100; %1020;         % density - Coating (kg/m^3)
mat.c_c       = 1900; %2480;         % speed of sound - Coating (m/s)
mat.E_c       = NaN;          % modulus of elasticity - Coating (Pa)
mat.sig_ult_c = 37e6;         % Ultimate Strength - Coating (Pa)
mat.m_c       = 16.92;%6.1;          % Wohler exponent - Coating
mat.nu_c      = 0.3;%0.42;         % Poisson's ratio - Coating

%Liquid - Water
mat.rho_l     = 1000;         % density - Liquid (kg/m^3)
mat.c_l       = 1480;         % speed of sound - Liquid (m/s)

%Material Thicknesses
mat.CoatTh    = 0.0004572; % Coating Thickness (m)
mat.LaminateTh= 0.0016;     % Leading Edge Laminate Thickness at tip (m)


[t_years_c, alpha_ltc, t_years_u, alpha_ltu, D_ltc, D_ltu] = ErosionDamageModelCoatedLaminate(phi_d, I, U_hub, P_I, stat, mat,...
                                                TurbineModel, RotorDiameter);
%%                                
[PercentMassLossCoating, PercentMassLossLaminate] = LEERollingHorizonSimCoatedLaminate(t_years_c, alpha_ltc,...
                                                    t_years_u, alpha_ltu, mat, BladeLifetime);

Time = linspace(0, BladeLifetime, BladeLifetime*365);
figure(1000); clf(1000)
plot(Time, 100*PercentMassLossCoating, 'LineWidth', 2)
hold on
plot(Time, 100*(PercentMassLossLaminate+PercentMassLossCoating), '--','LineWidth',1)
yline(100*0.1, '-', ' Category 1')
yline(100*0.5, '-', ' Category 2')
yline(  100*1, '-', ' Category 3')
yline(100*1.8, '-', ' Category 4')
yline(100*2.0, '-', ' Category 5')
xlabel('Turbine Lifetime (Years)')
ylabel('Mass Loss (%)')
ylim([0 210])
title('Leading edge Erosion (Inland)')
legend('Coating', 'Laminate', 'Location', 'northwest')


