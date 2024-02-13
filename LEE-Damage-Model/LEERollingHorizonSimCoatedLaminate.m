function [PercentMassLossCoating, PercentMassLossLaminate] = LEERollingHorizonSimCoatedLaminate(t_years_c, alpha_ltc, t_years_u, alpha_ltu, mat, BladeLifetime)
    
%     % Convert alpha_lt from kg/(m^2*s) ---> kg/(m^2*year)
%     alpha_lt_year = alpha_lt*31556926;
    
    % Convert alpha_lt from kg/(m^2*s) ---> kg/(m^2*day)
    alpha_ltc_day = alpha_ltc*86400;
    alpha_ltu_day = alpha_ltu*86400;
    
    %% Total normalized mass (kg/m^2)
    % Coating
    m_totalc = mat.rho_c*mat.CoatTh;
    % Laminate
    m_totalu = mat.rho_s*mat.LaminateTh;
    
    
    % Rolling Horizon simulation
    SimLength = BladeLifetime*365;
    FailureInitiationCoating = t_years_c*365;
    FailureInitiationLaminate = t_years_u*365;
    PercentMassLossCoating = zeros(1, length(SimLength)); %Initialize PercentMassLoss array
    PercentMassLossLaminate = zeros(1, length(SimLength)); %Initialize PercentMassLoss array
    CoatingDegradeTime = 0;
    
    for i = 1:SimLength % timestep = 1 day
        if i < FailureInitiationCoating
            PercentMassLossCoating(i) = 0;
        elseif PercentMassLossCoating(i-1) < 1
            m_c = alpha_ltc_day*(i - FailureInitiationCoating);
            PercentMassLossCoating(i) = m_c/m_totalc;
            CoatingDegradeTime = CoatingDegradeTime + 1;
        else 
            PercentMassLossCoating(i) = 1;
            if i > FailureInitiationLaminate
                % if time for laminate to fail exceeds time for coating to
                % fail + full coating degradation then use laminate failure
                % time to calculate mass loss, if it's less than coating
                % initiation + full failure use the summation to avoid
                % large discontinuity in mass loss prediction
                if FailureInitiationLaminate > (FailureInitiationCoating + CoatingDegradeTime) 
                    m_u = alpha_ltu_day*(i - FailureInitiationLaminate);
                    PercentMassLossLaminate(i) = m_u/m_totalu;
                    if PercentMassLossLaminate(i) >= 1
                        PercentMassLossLaminate(i) = 1;
                    end
                else
                    m_u = alpha_ltu_day*(i - (FailureInitiationCoating + CoatingDegradeTime));
                    PercentMassLossLaminate(i) = m_u/m_totalu;
                    if PercentMassLossLaminate(i) >= 1
                        PercentMassLossLaminate(i) = 1;
                    end
                end
                
            end
        end
            
    end
    
end