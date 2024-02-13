function [D_st, alpha_st] = ShortTermDamageRateSubstrate(phi_d, I, V_blade, mat)
    
    phi_d_m = phi_d/1000;
    % Calculate Erosive strength of material
    S = (4*mat.sig_ult_s*(mat.m_s - 1))/(1 - 2*mat.nu_sm); %Pa
    
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
    P_wh = (mat.rho_l*mat.c_l*V_imp)/(1 + (mat.rho_l*mat.c_l)/(mat.rho_s*mat.c_s)); % Pa
    
    % # of impacts till fatigue failure Ni
    ni = 7.1e-6*(S./P_wh).^5.7;
    
    % Dimensionless mass loss rate 
    alpha_star = 0.023*(1./ni).^0.7;
    
    % rate of mass loss (kg/impact)
    alpha_st = alpha_star.*(pi*mat.rho_s.*phi_d_m.^3)/4;
    
    % Multiply by impact rate to get temporal rate of mass loss (kg/m^2/s)
    alpha_st = alpha_st.*(q.*V_imp.*beta_d);
    
    % Short Term Damage Rate
    D_st = (q.*V_imp.*beta_d)./((8.9./phi_d.^2).*(S./P_wh).^5.7); % 1/s 
    
    
end