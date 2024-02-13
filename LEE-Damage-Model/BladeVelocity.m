function [TurbineModel] = BladeVelocity(TurbineModel, RotorDiameter)
    %% This function creates an array of blade velocities for each wind speed %%
    TipSpeed = [];
    C = pi*RotorDiameter;
    Wind = TurbineModel(:,1);
    for W = 1:length(Wind)
        RPM = TurbineModel(W,2);
        TipSpeed(W) = (C*RPM)/60; %meters/second
    end
    TurbineModel(:,3) = TipSpeed;
end