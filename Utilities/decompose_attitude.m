function [p_ia_dot, C_ia, C_ia_dot] = decompose_attitude(p_ia, w_ia_i)
%decompose_attitude Compute the DCM and the time rate of change of the DCM
%from vectorized DCM and the current angular velocity in the body frame.
%   Detail explanation goes here
%
% Ref: Keegan Bunker, Ryan Caverly, Modular Dynamic Modeling and Simulation
% of a Cable-Actuated Flexible Solar Sail, AIAA SCITECH 2024 Forum
%
% ----------------------------------------------------------------------- %

% Compute DCM
C_ia = p2DCM(p_ia);

% Compute d/dt DCM
C_ia_dot = -crossr(w_ia_i)*C_ia;

% Vectorize d/dt DCM
c1_dot = C_ia_dot(1,:); % Row 1
c2_dot = C_ia_dot(2,:); % Row 2
c3_dot = C_ia_dot(3,:); % Row 3

% Assemble vector
p_ia_dot = [c1_dot , c2_dot, c3_dot]'; % Row^T = Column

end