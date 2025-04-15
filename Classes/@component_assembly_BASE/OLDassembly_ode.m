function dydt = assembly_ode(obj, t, y)
%assembly_ode Computes the derivative for the general system based on the
%general null-space ODE.
%   Detailed explanation goes here
%
% TODO: This function has several components that are still hardcoded for
% the specific CABLESSail implementation used for SciTech. This must be
% more generalized, specifically the generalizes forces and nonlinear team
% need to be made generalized and then defined by a user defined file so
% this function can be used for any case

% ----------------------------------------------------------------------- %
    % Extract states
    n_q = sum(obj.state.dimensions);
    n_nu = sum(obj.velocity.dimensions);
    q  = y(1:n_q, end);
    nu_hat = y(n_q+1:n_q+n_nu,end);
    
    % Cell array of current states
    current_states = obj.state.decompose(q);

    % Cell array of current velocities
    current_velocities = obj.velocity.decompose(nu_hat);
    
    % Compute current null space matrix
    Upsilon     = obj.compute_upsilon(current_states, current_velocities);
    Upsilon_dot = obj.compute_upsilon_dot(current_states, current_velocities);
    
    % Compute pre-null-space states
    nu = Upsilon * nu_hat;
    
    % THIS WILL BREAK RIGHT HERE UNTIL I FIX UPSILON
    y_cell = obj.create_y_vectors(q, nu, nu_hat);

    % Current state matrices
    K       = obj.compute_K(t, y_cell);
    M_bar   = obj.compute_M_bar(t, y_cell);
    f_non   = obj.compute_f_non(t, y_cell);
    f       = obj.compute_f(t, y_cell);
    gamma   = obj.compute_Gamma(t, y_cell);

    % Compute natural frequencies
    % MinvK = inv(M_bar(13:18, 13:18))*K(25:30, 25:30);
    % wn2 = eig(MinvK);
    % wnq = 1/(2*pi)*sqrt(wn2);
    
    %-------------------------------------------------------------%
    % REACTION FORCES
    %-------------------------------------------------------------%
    % Set generalized force vectors to 0
    f_r_ba_a = [0 0 0]';
    f_w_ba_b = [0 0 0]';
    f_q1 = zeros(6,1);
    f_q2 = zeros(6,1);
    f_q3 = zeros(6,1);
    f_q4 = zeros(6,1);

    % DCM's between booms and spacecraft bus
    C_il = @(qii) [1 0 -obj.components{2}.const.dpsi_i_dx_FH(29.5)*qii;...
               0 1 0;...
               -obj.components{2}.const.dpsi_i_dx_FH(29.5)*qii 0 1];
    C_cb = eye(3);
    
    C_db = [0 1 0;
           -1 0 0;
           0 0 1];
    C_eb = [-1 0 0;
            0 -1 0;
            0 0 1];
    C_fb = [0 -1 0;
            1 0 0;
            0 0 1];
    
    % Position of boom bases relative to bus cg
    R = 1; % meter
    r_cb_b = R*[1;0;0];
    r_db_b = R*[0;1;0];
    r_eb_b = R*[-1;0;0];
    r_fb_b = R*[0;-1;0];

    % TODO put the reaction forces in the component level function
    % Extract elastic coordinates to compute generalized forces
    q11 = current_states{5}(1:3);
    q12 = current_states{5}(4:6);
    q1 = [q11;q12];

    q21 = current_states{8}(1:3);
    q22 = current_states{8}(4:6);
    q2 = [q21;q22];

    q31 = current_states{11}(1:3);
    q32 = current_states{11}(4:6);
    q3 = [q31;q32];

    q41 = current_states{14}(1:3);
    q42 = current_states{14}(4:6);
    q4 = [q41;q42];


    % Cable Reaction forces --------------------------------------
    %-------------------------------------------------------------%
    
    % X values where plate locations are
    % Compute cable reaction forces
    F11 = -compute_cable_forces(T11, obj.components{2}.const.x_plates, q11, obj.components{2}.const);
    F12 = -compute_cable_forces(T12, obj.components{2}.const.x_plates, q12, obj.components{2}.const);

    F21 = -compute_cable_forces(T21, obj.components{3}.const.x_plates, q21, obj.components{3}.const);
    F22 = -compute_cable_forces(T22, obj.components{3}.const.x_plates, q22, obj.components{3}.const);

    F31 = -compute_cable_forces(T31, obj.components{4}.const.x_plates, q31, obj.components{4}.const);
    F32 = -compute_cable_forces(T32, obj.components{4}.const.x_plates, q32, obj.components{4}.const);

    F41 = -compute_cable_forces(T41, obj.components{5}.const.x_plates, q41, obj.components{5}.const);
    F42 = -compute_cable_forces(T42, obj.components{5}.const.x_plates, q42, obj.components{5}.const);

    % Effect on elastic coordinates
    f_q1_cable_reaction = zeros(6,1);
    f_q2_cable_reaction = zeros(6,1);
    f_q3_cable_reaction = zeros(6,1);
    f_q4_cable_reaction = zeros(6,1);
    for lv1 = 1:length(obj.components{2}.const.x_plates)-2
        xi = obj.components{2}.const.x_plates(lv1+1);
        f_q1_cable_reaction = f_q1_cable_reaction + ([0; F11(lv1); F12(lv1)]'*C_cb'*[zeros(1,6); obj.components{2}.const.PSI_FH(xi)])';
        f_q2_cable_reaction = f_q2_cable_reaction + ([0; F21(lv1); F22(lv1)]'*C_cb'*[zeros(1,6); obj.components{3}.const.PSI_FH(xi)])';
        f_q3_cable_reaction = f_q3_cable_reaction + ([0; F31(lv1); F32(lv1)]'*C_cb'*[zeros(1,6); obj.components{4}.const.PSI_FH(xi)])';
        f_q4_cable_reaction = f_q4_cable_reaction + ([0; F41(lv1); F42(lv1)]'*C_cb'*[zeros(1,6); obj.components{5}.const.PSI_FH(xi)])';
    end

    F_cable_reaction = [f_q1_cable_reaction;
                        f_q2_cable_reaction;
                        f_q3_cable_reaction;
                        f_q4_cable_reaction];

    % Ignore the effect on omega and position since those generalized
    % forces are canceled since these are internal forces
    %-------------------------------------------------------------%

    % Cable Reaction tip force -----------------------------------
    %-------------------------------------------------------------%
    % TODO: Move this to the boom component as a method instead. These
    % computations only depend on the current tension and current elastic
    % coordinate state. 
    
    % f_q12_cable_momet = ( (C_il(q12)*[-T12; 0; 0])'*C_cb'*[-(obj.components{2}.const.Radius + obj.components{2}.const.Plate_h)*obj.components{2}.const.dpsi_i_dx_FH(29.5); zeros(1,3); obj.components{2}.const.psi_i_FH(29.5)] )';
    % f_q22_cable_momet = ( (C_il(q22)*[-T22; 0; 0])'*C_db'*[-(obj.components{3}.const.Radius + obj.components{3}.const.Plate_h)*obj.components{3}.const.dpsi_i_dx_FH(29.5); zeros(1,3); obj.components{3}.const.psi_i_FH(29.5)] )';
    % f_q32_cable_momet = ( (C_il(q32)*[-T32; 0; 0])'*C_eb'*[-(obj.components{4}.const.Radius + obj.components{4}.const.Plate_h)*obj.components{4}.const.dpsi_i_dx_FH(29.5); zeros(1,3); obj.components{4}.const.psi_i_FH(29.5)] )';
    % f_q42_cable_momet = ( (C_il(q42)*[-T42; 0; 0])'*C_fb'*[-(obj.components{5}.const.Radius + obj.components{5}.const.Plate_h)*obj.components{5}.const.dpsi_i_dx_FH(29.5); zeros(1,3); obj.components{5}.const.psi_i_FH(29.5)] )';

    % TODO: These need to be added for the inplane bending booms as well
    f_q12_cable_momet = T12 * (obj.components{2}.const.Radius + obj.components{2}.const.Plate_h) * obj.components{2}.const.dpsi_i_dx_FH(29.5);
    f_q22_cable_momet = T22 * (obj.components{3}.const.Radius + obj.components{3}.const.Plate_h) * obj.components{3}.const.dpsi_i_dx_FH(29.5);
    f_q32_cable_momet = T32 * (obj.components{4}.const.Radius + obj.components{4}.const.Plate_h) * obj.components{4}.const.dpsi_i_dx_FH(29.5);
    f_q42_cable_momet = T42 * (obj.components{5}.const.Radius + obj.components{5}.const.Plate_h) * obj.components{5}.const.dpsi_i_dx_FH(29.5);
    
    F_cable_moment = [zeros(3,1);
                     f_q12_cable_momet';
                     zeros(3,1);
                     f_q22_cable_momet';
                     zeros(3,1);
                     f_q32_cable_momet';
                     zeros(3,1);
                     f_q42_cable_momet'];


    % Compute  SRP forces and torques ----------------------------
    %-------------------------------------------------------------%
    % Bus cg location wrt inertial poit a
    r_ba_a = current_states{1};

    % Boom tip deflections (point t) in the boom frames
    r_tc_c = [(obj.components{2}.const.Length); obj.components{2}.const.PSI_FH(obj.components{2}.const.Length)*q1];
    r_td_d = [(obj.components{3}.const.Length); obj.components{3}.const.PSI_FH(obj.components{3}.const.Length)*q2];
    r_te_e = [(obj.components{4}.const.Length); obj.components{4}.const.PSI_FH(obj.components{4}.const.Length)*q3];
    r_tf_f = [(obj.components{5}.const.Length); obj.components{5}.const.PSI_FH(obj.components{5}.const.Length)*q4];

    % Tip positions in Bus frame (sail body frame)
    r_tc_b = r_cb_b + C_cb' * r_tc_c;
    r_td_b = r_db_b + C_db' * r_td_d;
    r_te_b = r_eb_b + C_eb' * r_te_e;
    r_tf_b = r_fb_b + C_fb' * r_tf_f;

    % Sun direction in the inertial frame
    
    % Bus location
    Pos_val.hub = [0 0 0]';

    % Optical properties of the sail (Solar Cruiser parameters from Heaton and
    % Artusio-Glimpse "An Update to the NASA Reference Solal Sail Thrust
    % Model'' 2015
    Opt_prop.P = 4.53914e-6;
    Opt_prop.r_tilde = 0.91;
    Opt_prop.s       = 0.94;
    Opt_prop.Bf      = 0.79;
    Opt_prop.Bb	     = 0.67;
    Opt_prop.ef      = 0.025;
    Opt_prop.eb       = 0.27;

    % Form C_ab from time history
    p_ba_t = current_states{2};
    [C_ab, ~] = build_C_ab(p_ba_t);

    % Transform sun into the body frame
    % Assuming the be frame is aligne with the a grame at t=0
    % Use a 2-rotation to set the SIA in the plane of booms 1 and 3
    Sun_dir_sun_frame = [0;0;1];
    SIA = 0; % degrees
    % CHANGE IN SRP PLOT SCRIPT
    C_a_sun = [cosd(SIA) 0 sind(SIA);
              0 1 0;
              -sind(SIA) 0 cosd(SIA)];
    % assuming C_ab = eye(3) at t=0
    Sun_dir_b = C_ab' * C_a_sun*Sun_dir_sun_frame;
    
    % Compute force and torque from SRP

    % Format for sail code
    % Pos_val.tip_c = r_tc_b;
    % Pos_val.tip_d = r_td_b;
    % Pos_val.tip_e = r_te_b;
    % Pos_val.tip_f = r_tf_b;

    % Sail 1
    Pos_val.tip1 = r_tc_b;
    Pos_val.tip2 = r_td_b;
    [Force_sail1,Torque_sail1] = RigidSailModel(Pos_val,Sun_dir_b,Opt_prop);
    
    % Sail 2
    Pos_val.tip1 = r_td_b;
    Pos_val.tip2 = r_te_b;
    [Force_sail2,Torque_sail2] = RigidSailModel(Pos_val,Sun_dir_b,Opt_prop);
    
    % Sail 3
    Pos_val.tip1 = r_te_b;
    Pos_val.tip2 = r_tf_b;
    [Force_sail3,Torque_sail3] = RigidSailModel(Pos_val,Sun_dir_b,Opt_prop);
    
    % Sail 4
    Pos_val.tip1 = r_tf_b;
    Pos_val.tip2 = r_tc_b;
    [Force_sail4,Torque_sail4] = RigidSailModel(Pos_val,Sun_dir_b,Opt_prop);
    
    % Total force and torques
    Force_sail = Force_sail1.total + Force_sail2.total + Force_sail3.total + Force_sail4.total;
    Torque_sail = Torque_sail1.total + Torque_sail2.total + Torque_sail3.total + Torque_sail4.total;

    % External Force ---------------------------------------------
    %-------------------------------------------------------------%
    F_ext_a = Force_sail;

    % External Moment -------------------------------------------
    %-------------------------------------------------------------%
    T_ext = Torque_sail;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assembe total generalized force vector

    F_gen_r = F_ext_a;
    F_gen_w = T_ext;

    % Turn things off -----------------------
        F_gen_r = zeros(3,1);
        F_gen_w = zeros(3,1);
        %F_cable_moment = zeros(4*6,1);
        %F_cable_reaction = zeros(4*6,1);
    % ---------------------------------------
    
    F_gen = [F_gen_r; F_gen_w; F_cable_moment+F_cable_reaction];
    
    % if t>10
    %     disp('help')
    % end

    %F_gen = zeros(30,1);
    %F_gen(3) = -1;
    % Add damping ------------------------------------------------
    % 
    damp_ratio = 0.00001;
    % damp_ratio = 0;
    w_n_q = [3 1 3]';
    D_q = diag([damp_ratio*w_n_q; damp_ratio*w_n_q]);
    D = blkdiag(zeros(3,3), zeros(3,3), D_q, D_q, D_q, D_q);

    % ----------------------------------------
    % Checking the actual damping of the system
    a11 = zeros(84,84);
    a12 = gamma*Upsilon;
    a21 = -inv(Upsilon'* M_bar * Upsilon)*Upsilon'*gamma'*K;
    a22 = -inv(Upsilon'* M_bar * Upsilon)*D;
    A_linear = [a11, a12;
                a21, zeros(30,30)];
    A_linear_damp = [a11, a12;
                a21, a22];
    w_n2 = abs(eig(A_linear));
    zeta = -cos(angle(eig(A_linear_damp)));
    % ----------------------------------------



    % Update state
    v_dot = inv(Upsilon'* M_bar * Upsilon) * (F_gen - D*nu_hat+ Upsilon'*(gamma'*f - gamma'*K*q - f_non - M_bar*Upsilon_dot*nu_hat));
    dydt = [gamma*Upsilon*nu_hat; v_dot];
end