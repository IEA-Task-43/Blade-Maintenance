function [D_st, alpha_st] = ShortTermDamageRateCoatedLaminate(phi_d, I, V_blade, mat)
    
    %Convert phi_d from mm to meters
    phi_d_m = phi_d/1000;
    
    % Calculate impedances
    Z_s = mat.rho_s*mat.c_s; % Substrate impedance
    Z_l = mat.rho_l*mat.c_l; % Liquid impedance
    Z_c = mat.rho_c*mat.c_c; % Coating impedance
    
    % Calculate Relative Impedances
    psi_sc = (Z_s - Z_c)/(Z_s + Z_c);
    psi_lc = (Z_l - Z_c)/(Z_l + Z_c);
    
    % quantity representing the proportional thickness of the coating to
    % droplet with relative impedances
%     gamma = (mat.c_c + Z_l/Z_s)/(mat.c_l + Z_c/Z_s)*(2/(1+Z_l/Z_c))*(phi_d_m/mat.CoatTh);
    gamma = (mat.c_c/mat.c_l)*(phi_d_m/mat.CoatTh)*((1 + (Z_l/Z_s))/(1 + (Z_c/Z_s)))*((1 + Z_l/Z_s)/2);
    
    %# of wave reflections that occur while the droplet is impacting the coating
    k = (1 - exp(-gamma))/(1 - psi_sc*psi_lc);
                                
    % Calculate Erosive strength of coating material
    S_c = (4*mat.sig_ult_c*(mat.m_c - 1))/(1 - 2*mat.nu_c); %Pa  
    
    % Modified Strength of coating
    S_ec = S_c./(1 + 2*k'*abs(psi_sc));
    
    %% OG Code
    % Calculate terminal speed of individual rain droplets
    
    V_tg = 9.65 - 10.30*exp(-0.6*phi_d); % m/s

    
    % Calculate Impact velocity
    % Calculate number of rain droplets/unit volume of rainfall
    x = repmat(V_blade,length(V_tg),1,1);
    y = repmat(V_tg,1,1,length(V_blade));
    V_imp = x + y; % m/s
    q = 530.5.*(I./(V_tg.*phi_d.^3));

    V_imp(V_imp<0) = 0;
    q(q<0) = 0;
    
    % Calculate Impingement efficiency
    beta_d = 1 - exp(-15.*phi_d); % unitless
    
    % Calculate Water Hammer Pressure 
    P_wh = (mat.rho_l*mat.c_l*V_imp)/(1 + (mat.rho_l*mat.c_l)/(mat.rho_c*mat.c_c)); % Pa
    %%
    % New Pressure felt at coating interface
    sig_o = (1 + psi_sc)/(1 - psi_sc*psi_lc)*(1 - psi_sc*((1 + psi_lc)/(1 + psi_sc))*((1 - exp(-gamma))/gamma)).*P_wh;
    
    % # of impacts till fatigue failure Ni
    ni = 7.1e-6.*(S_ec./sig_o).^5.7;
    
    % Dimensionless mass loss rate 
    alpha_star = 0.023*(1./ni).^0.7;
    
    % rate of mass loss (kg/impact)
    alpha_st = alpha_star.*(pi*mat.rho_c.*phi_d_m.^3)/4;
    
    % Multiply by impact rate to get temporal rate of mass loss (kg/m^2/s)
    alpha_st = alpha_st.*(q.*V_imp.*beta_d);
    
    % Short Term Damage Rate
    D_st = (q.*V_imp.*beta_d)./((8.9./phi_d.^2).*(S_ec./sig_o).^5.7); % 1/s 
    
   c=1;
   
end