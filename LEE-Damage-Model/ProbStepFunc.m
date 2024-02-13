function PI_step = ProbStepFunc(P_I, I)
    %% Creating step function vector to represent changing Probabilities for Rain Intensity%%
    
    PI_step = zeros(1,length(I));
    PI_step(I>=P_I(1).Min & I<P_I(1).Max) = P_I(1).discretizedProbability;
    PI_step(I>=P_I(2).Min & I<P_I(2).Max) = P_I(2).discretizedProbability;
    PI_step(I>=P_I(3).Min & I<P_I(3).Max) = P_I(3).discretizedProbability;
    PI_step(I>=P_I(4).Min)                = P_I(4).discretizedProbability;
    
end