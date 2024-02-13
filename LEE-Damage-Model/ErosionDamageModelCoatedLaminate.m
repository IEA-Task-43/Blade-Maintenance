function [t_years_c, alpha_ltc, t_years_u, alpha_ltu, D_ltc, D_ltu] = ErosionDamageModelCoatedLaminate(phi_d, I, U_hub, P_I, stat, mat,...
                                                TurbineModel, RotorDiameter)
%% 
% I             = Range of Rainfall intensities
% phi_d         = Range of Rain droplet diameters
% U_hub         = Range of wind speeds at the hub of the turbine
% P_I           = Structure of rainfall condition ranges and their associated probabilities

% mu_I          = Mean of rainfall intensity
% sig_I         = Standard deviation of rainfall intensity
% alpha_u       = Weibull shape parameter for wind speed distribution
% beta_u        = Weibull scale parameter for wind speed distribution
% N             = Weibull parameter for Droplet Size Distribution
% q             = Weibull parameter for Droplet Size Distribution
% A             = Weibull parameter for Droplet Size Distribution
% p             = Weibull parameter for Droplet Size Distribution

%% Laminate
% rho_s         = density - Laminate (kg/m^3)
% c_s           = speed of sound- Laminate (m/s)
% E_sm          = Modulus of elasticity - Laminate matrix (Pa)
% sig_ult_s     = Laminate matrix material ultimate strength (Pa)
% m_s           = Wohler fatigue slope of laminate matrix material
% nu_sm         = Poisson's ratio of Laminate Matrix material
% sige_sm       = Endurance Limit Laminate Matrix material (Pa)
% b_s           = Fatigue Kneee Laminate Matrix material
% Ell_s         = Young's Modulus in fiber direction - Laminate (Pa)
% E22_s         = Young's Modulus in transverse direction - Laminate (Pa)
% G12_s         = Shear Modulus - Laminate (Pa)
% nu12_s        = Poisson's ratio - Laminate
% LaminateTh    = Laminate Thickness at Leading edge

%% Coating
% rho_c         = Density - Coating (kg/m^3)
% c_c           = Speed of Sound - Coating (m/s)
% E_c           = Modulus of Elasticity - Coating (Pa)
% sig_ult_c     = Ultimate Strength - Coating (Pa)
% m_c           = Wohler exponent - Coating
% nu_c          = Poisson's ratio - Coating
% CoatTh        = Coating Thickness (m)

%% Liquid - Water
% rho_w         = density of water, set at 1000kg/m^3
% c_w           = speed of sound of water, set at  1480 m/s 

% TurbineModel  = Pandas dataframe of wind speed vs. blade velocity at some span, usually the tip
% RotorDiameter = Length at which LEE is being analyzed
%%

alpha_u = stat.alpha_u;
beta_u  = stat.beta_u;
mu_I    = stat.mu_I;
sig_I   = stat.sig_I;
N       = stat.N; 
q       = stat.q;
A       = stat.A; 
p       = stat.p;

% sig_ult = mat.sig_ult; % Pa;
% m       = mat.m;
% nu      = mat.nu;
% rho_s   = mat.rho_s;   % kg/m^3
% c_s     = mat.c_s;     % m/s
% rho_w   = mat.rho_w;   % kg/m^3
% c_w     = mat.c_w;     % m/s

TurbineModel = BladeVelocity(TurbineModel, RotorDiameter);
Cut_in = TurbineModel(1,1);     %Determine cut-in wind speed
Cut_out = TurbineModel(end,1);  %Determine cut-out wind speed

%Determine # of discretizations contained within each rainfall condition
%range
for x = 1:length(P_I)
    if isnan(P_I(x).Max)
        %Set Max value to very large number
        P_I(x).Max = 1000000;
    end
    P_I(x).discretizedProbability = (P_I(x).Probability/100);%/counter;
end

% Make probability vector
PI_step = ProbStepFunc(P_I, I);

%Interpolate blade velocity to discretized wind speed
U_tip = zeros(1, length(U_hub));

for x = 1:length(U_hub)
    if U_hub(x) < Cut_in 
        U_tip(x) = 0;
    elseif U_hub(x) > Cut_out
        U_tip(x) = 0;
    else
        U_tip(x) = interp1(TurbineModel(:,1), TurbineModel(:,3), U_hub(x));
    end
end

% Reshape vectors for multiplication %%%%%%%%%%%%%%%
U_hub2 = reshape(U_hub, [1,1,length(U_hub)]);
U_tip = reshape(U_tip, [1,1,length(U_tip)]);
phi_d = phi_d';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSD_PDF = createDSD_General(phi_d, I, N, q, A, p);

I_PDF = createRainIntensityDistribution(I, mu_I, sig_I);
            
f_I_phi = I_PDF.*DSD_PDF; % Calculate Joint PDF of droplet size and rainfall intensity
           
U_PDF = createWindDistribution(U_hub2, alpha_u, beta_u); % Wind statistics model
           
% short term damage rates and mass loss rates 
[D_stc, alpha_stc] = ShortTermDamageRateCoatedLaminate(phi_d, I, U_tip, mat);
[D_stu, alpha_stu] = ShortTermDamageRateSubstrate(phi_d, I, U_tip, mat);
                                      

% Coated Laminate
DLTCfunc = D_stc.*f_I_phi.*PI_step.*U_PDF;
AlphaLTCfunc = alpha_stc.*f_I_phi.*PI_step.*U_PDF;
D_ltc = trapz(phi_d, trapz(I, trapz(U_hub, DLTCfunc, 3), 2), 1);
%Alpha_lt is in kg/m^2/s
alpha_ltc = trapz(phi_d, trapz(I, trapz(U_hub, AlphaLTCfunc, 3), 2), 1);


% Uncoated Laminate
DLTUfunc = D_stu.*f_I_phi.*PI_step.*U_PDF;
AlphaLTUfunc = alpha_stu.*f_I_phi.*PI_step.*U_PDF;
D_ltu = trapz(phi_d, trapz(I, trapz(U_hub, DLTUfunc, 3), 2), 1);
%Alpha_lt is in kg/m^2/s
alpha_ltu = trapz(phi_d, trapz(I, trapz(U_hub, AlphaLTUfunc, 3), 2), 1);
                    
%% Convert long-term damage rate to coating initiation lifetime
%Coated Laminate
t_years_c = 1/(D_ltc*(365*24*3600));
%Uncoated Laminate
t_years_u = 1/(D_ltu*(365*24*3600));

end


