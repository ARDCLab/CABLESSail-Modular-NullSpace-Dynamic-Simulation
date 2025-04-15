function f_non = compute_f_non(obj, t, coords, vels, q, v)
    %compute_f_non Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ----------------------------------------------------------- %
    % Extract states and velocities
    velocity = [vels{1}; vels{2}];
    r_dot = vels{1};
    omega = vels{2};
    C_ia = p2DCM(coords{2});
    %[~, C_ia, ~] = decompose_attitude(coords{2}, omega);
    
    % Assemble component matrices and vectors
    %Omega_x = obj.compute_Omega(omega);
    Omega_x = blkdiag( zeros(length(velocity)-3), crossr(omega));
    M       = obj.compute_M(t, coords, vels, q, v);
    M_dot   = M*0; % Mdot = 0  for the hub?
    
    % Compute f_nonlinear
    f_non = M_dot*velocity + Omega_x * M * velocity;


    %f_non = zeros(6,1);
end