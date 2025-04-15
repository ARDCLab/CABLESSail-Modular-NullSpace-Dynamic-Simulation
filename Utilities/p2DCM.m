function [C_ia] = p2DCM(p_ia)
%p2DCM Convert the vectorized DCM state into the DCM
%   Detailed explanation goes here
%
% Ref: Keegan Bunker, Ryan Caverly, Modular Dynamic Modeling and Simulation
% of a Cable-Actuated Flexible Solar Sail, AIAA SCITECH 2024 Forum
%
% ----------------------------------------------------------------------- %
    p1 = p_ia(1:3);
    p2 = p_ia(4:6);
    p3 = p_ia(7:9);
    C_ai = [p1, p2, p3];
    C_ia = C_ai';
end