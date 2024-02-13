function I_PDF = createRainIntensityDistribution(I, mu_I, sig_I)
    %% Assumed Lognormal distribution for rainfall intensity PDF %%
    
    I_PDF = 1./(sqrt(2.*pi).*I.*sig_I).*exp(-(log(I) - mu_I).^2/(2.*sig_I.^2));
    
end