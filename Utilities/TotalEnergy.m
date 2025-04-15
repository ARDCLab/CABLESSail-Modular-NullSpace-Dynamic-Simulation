function E = TotalEnergy(t,q, q_dot, M, K)
% Calculate potential, kinetic, and change in total energy over time.
%
% Inputs
%   t -     Vector of time steps. (1 x n)
%   q -     Matrix of all states through time. (m x n) m = # of states
%   q_dot - Matrix of state time derivativs. (m x n)
%   M -     Mass matrix
%   K -     Stiffness matrix
%
% Outputs
%   E -     Matrix of all energies through time. (3 x n)
%           E(1,:) - kinetic energy
%           E(2,:) - potential energy
%           E(3,:) - change in total energy
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: May 31, 2023
%-------------------------------------------------------------------------%
% TODO:
%   -   Calculate total energy as well, instead of just the change?
%
%
%-------------------------------------------------------------------------%

E = zeros(3,max(size(q)));
E0 = 0;
for lv1 = 1:max(size(q))
    qi = q(:,lv1);
    qi_dot = q_dot(:,lv1);
    % Kinetic energy
    E(1,lv1) = 1/2 * qi_dot'*M*qi_dot;
    % Potential energy
    E(2,lv1) = 1/2 * qi'*K*qi;
    % Total energy
    E(3,lv1) = 1/2 * qi_dot'*M*qi_dot + 1/2 * qi'*K*qi; %sum(E(:,lv1));
    if lv1==1
        E0 = E(3,lv1);
    end
    E(3,lv1) = E(3,lv1)- E0;
end
end