function ups_dot = compute_Upsilon_dot(obj, t, coords, vels_hat, q, v_hat)  
%compute_upsilon_dot User-defined function that computes the total Upsilon
%matrix of the assembled system.
%
%   Detailed explanation goes here
%
% ----------------------------------------------------------------------- %
    % Form DCM
    %[C_ab, ~] = build_C_ab(coords{2}, false);
    C_ba = p2DCM(coords{2});
    %[C_ia] = p2DCM(p_ia)
    C_ab = C_ba';

    % Extract omega
    % TODO: check w_ba/w_ab
    % From the reference, the state is w^ba_b, so this should be fine
    w_ba_b = vels_hat{2};

    % Get the boom attachment locations
    r_cb_b = obj.const.r_cb_b;
    r_db_b = obj.const.r_db_b;
    r_eb_b = obj.const.r_eb_b;
    r_fb_b = obj.const.r_fb_b;
    
    % Compute Upsilon_dot
    ups_dot_c = -C_ab*crossr(w_ba_b)*crossr(r_cb_b);
    ups_dot_d = -C_ab*crossr(w_ba_b)*crossr(r_db_b);
    ups_dot_e = -C_ab*crossr(w_ba_b)*crossr(r_eb_b);
    ups_dot_f = -C_ab*crossr(w_ba_b)*crossr(r_fb_b);

    ups_dot = [zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    ups_dot_c,  zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3), zeros(6), zeros(6), zeros(6), zeros(6);
               zeros(3),    ups_dot_d,  zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3), zeros(6), zeros(6), zeros(6), zeros(6);
               zeros(3),    ups_dot_e,  zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3), zeros(6), zeros(6), zeros(6), zeros(6);
               zeros(3),    ups_dot_f,  zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    zeros(3),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3), zeros(6), zeros(6), zeros(6), zeros(6)];
end
