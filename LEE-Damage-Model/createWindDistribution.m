function U_PDF = createWindDistribution(U_hub, alpha_u, beta_u)
    %% Assumed two-parameter weibull distribution for rainfall intensity PDF %%
    
    U_PDF = (alpha_u./beta_u).*(U_hub./beta_u).^(alpha_u - 1).*exp(-(U_hub./beta_u).^alpha_u);
    
end