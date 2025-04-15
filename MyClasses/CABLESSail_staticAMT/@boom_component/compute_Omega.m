% Omega matrix function handle, input is angular velocity
function Omega_x = compute_Omega(obj, omega)
%compute_Omega Summary of this function goes here
%   Detailed explanation goes here
%
% ------------------------------------------------------------------------%
    Omega_x = blkdiag(zeros(3,3), crossr(omega), zeros(obj.velocity.states{3}.dimension) );
end