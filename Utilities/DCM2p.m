function [p_ia] = DCM2p(C_ia)
%p2DCM Convert the vectorized DCM state into the DCM
%   Detailed explanation goes here
%
% Ref: Keegan Bunker, Ryan Caverly, Modular Dynamic Modeling and Simulation
% of a Cable-Actuated Flexible Solar Sail, AIAA SCITECH 2024 Forum
%
% ----------------------------------------------------------------------- %
 
C_ai = C_ia';
p1 = C_ai(:,1);
p2 = C_ai(:,2);
p3 = C_ai(:,3);
p_ia = [p1; p2; p3];

end