function Upsilon = compute_Upsilon(obj, t, coords, vels_hat, q, v_hat)  
%compute_upsilon User-defined function that computes the total Upsilon
%matrix of the assembled system.
%
%   Detailed explanation goes here

% ----------------------------------------------------------------------- %

    % Form DCM
    %[C_ab, ~] = build_C_ab(coords{2}, false);
    C_ba = p2DCM(coords{2});
    %[C_ia] = p2DCM(p_ia)
    C_ab = C_ba';

    % Get the boom orientations
    C_cb = obj.const.C_cb;
    C_db = obj.const.C_db;
    C_eb = obj.const.C_eb;
    C_fb = obj.const.C_fb;

    % Get the boom attachment locations
    r_cb_b = obj.const.r_cb_b;
    r_db_b = obj.const.r_db_b;
    r_eb_b = obj.const.r_eb_b;
    r_fb_b = obj.const.r_fb_b;

    % Get the AMT position
    r_gb_b = obj.const.r_gb_b;

    % Compute Upsilon
    Upsilon = [eye(3),      zeros(3),               zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    eye(3),                 zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);

               eye(3),      -C_ab*crossr(r_cb_b),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    C_cb,                   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3),             eye(6),     zeros(6),   zeros(6),   zeros(6);

               eye(3),      -C_ab*crossr(r_db_b),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    C_db,                   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3),             zeros(6),   eye(6),     zeros(6),   zeros(6);

               eye(3),      -C_ab*crossr(r_eb_b),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    C_eb,                   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3),             zeros(6),   zeros(6),   eye(6),     zeros(6);

               eye(3),      -C_ab*crossr(r_fb_b),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    C_fb,                   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(6,3),  zeros(6,3),             zeros(6),   zeros(6),   zeros(6),   eye(6);
               
               eye(3),      -C_ab*crossr(r_gb_b),   zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6);
               zeros(3),    eye(3),                 zeros(3,6), zeros(3,6), zeros(3,6), zeros(3,6)];

end
