function f_non = compute_f_non(obj, t, coords, vels, q, v)
    %compute_f_non Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ----------------------------------------------------------- %
    % Extract states and velocities
    velocity = [vels{1}; vels{2}; vels{3}];
    r_dot = vels{1};
    omega = vels{2};
    q_eps_dot = vels{3};
    C_ia = p2DCM(coords{2});
    %[~, C_ia, ~] = decompose_attitude(coords{2}, omega);
    
    % Assemble component matrices and vectors
    Omega_x = obj.compute_Omega(omega);
    M_dot   = obj.compute_Mdot(t, coords, vels, q, v);
    M       = obj.compute_M(t, coords, vels, q, v);
    v_bi    = obj.const.H_B_i*q_eps_dot - obj.const.c_B_ix * omega;
    
    % Compute f_nonlinear
    f_non = M_dot*velocity + Omega_x * M * velocity - [zeros(3,1); crossr(v_bi)*C_ia*r_dot; zeros(6,1)];
end 