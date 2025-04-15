function obj = set_initial_conditions(obj,coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% ----------------------------------------------------------------------- %

% ------------------------------------------------------------ %
%   TODO?: Make a function to go from q_reduced --> q          %
% ------------------------------------------------------------ %
% 
% % Get DCM
% C_ba = p2DCM(p_ba);
% C_ab = C_ba';
% 
% % Get the boom orientations
% C_cb = obj.const.C_cb;
% C_db = obj.const.C_db;
% C_eb = obj.const.C_eb;
% C_fb = obj.const.C_fb;
% 
% % Get the boom attachment locations
% r_cb_b = obj.const.r_cb_b;
% r_db_b = obj.const.r_db_b;
% r_eb_b = obj.const.r_eb_b;
% r_fb_b = obj.const.r_fb_b;
% 
% % Resolve boom locations in the inerital frame
% r_cb_a = C_ab* r_cb_b;
% r_db_a = C_ab* r_db_b;
% r_eb_a = C_ab* r_eb_b;
% r_fb_a = C_ab* r_fb_b;
% 
% r_ca_a = r_cb_a  + r_ba_a;
% r_da_a = r_db_a  + r_ba_a;
% r_ea_a = r_eb_a  + r_ba_a;
% r_fa_a = r_fb_a  + r_ba_a;
% 
% % Compute DCMs of booms
% C_ca = C_cb*C_ba;
% C_da = C_db*C_ba;
% C_ea = C_eb*C_ba;
% C_fa = C_fb*C_ba;
% 
% % Decompose DCM's
% [p_ca] = DCM2p(C_ca);
% [p_da] = DCM2p(C_da);
% [p_ea] = DCM2p(C_ea);
% [p_fa] = DCM2p(C_fa);
% 
% % Form state vectors
% q_Bb_0 = [r_ba_a; p_ba];
% q_Bc_0 = [r_ca_a; p_ca; q1];
% q_Bd_0 = [r_da_a; p_da; q2];
% q_Be_0 = [r_ea_a; p_ea; q3];
% q_Bf_0 = [r_fa_a; p_fa; q4];
% 
% % Set initial condtions of state vectors in component objects
% obj.components{1}.coordinate.set_initial_condition(q_Bb_0);
% obj.components{2}.coordinate.set_initial_condition(q_Bc_0);
% obj.components{3}.coordinate.set_initial_condition(q_Bd_0);
% obj.components{4}.coordinate.set_initial_condition(q_Be_0);
% obj.components{5}.coordinate.set_initial_condition(q_Bf_0);
% 
% % Compute all states and velocities
% q_0 = [q_Bb_0; q_Bc_0; q_Bd_0; q_Be_0; q_Bf_0];
% v_hat_0 = [v_ba; w_ba; q1dot; q2dot; q3dot; q4dot];
% y_0 = [q_0; v_hat_0];
% ------------------------------------------------------------ %


% Compute the expanded set of coordiantes
[coords_0, q_0] = obj.q_hat2q(coords_reduced_0, q_hat_0);

% Things are redefined a lot, put in breakpoints and check coordinate
% values to make sure they stay consistent


% Compute the expanded set of velocities
y_0 = [q_0; v_hat_0];
[coords_0, vels_0, vels_reduced_0, q_0, v_0, v_hat_0, Upsilon, Upsilon_dot] = obj.decomepose_state(y_0);
[comp_coords_array_0, comp_vels_array_0, comp_q_array, comp_v_array] = decompoese_component_states(obj, coords_0, vels_0, q_0, v_0);

% This should give the same answer
obj.components{1}.coordinate = obj.components{1}.coordinate.set_initial_condition(comp_q_array{1});
obj.components{2}.coordinate = obj.components{2}.coordinate.set_initial_condition(comp_q_array{2});
obj.components{3}.coordinate = obj.components{3}.coordinate.set_initial_condition(comp_q_array{3});
obj.components{4}.coordinate = obj.components{4}.coordinate.set_initial_condition(comp_q_array{4});
obj.components{5}.coordinate = obj.components{5}.coordinate.set_initial_condition(comp_q_array{5});

% Set initial conditions of component velocities
obj.components{1}.velocity = obj.components{1}.velocity.set_initial_condition(comp_v_array{1});
obj.components{2}.velocity = obj.components{2}.velocity.set_initial_condition(comp_v_array{2});
obj.components{3}.velocity = obj.components{3}.velocity.set_initial_condition(comp_v_array{3});
obj.components{4}.velocity = obj.components{4}.velocity.set_initial_condition(comp_v_array{4});
obj.components{5}.velocity = obj.components{5}.velocity.set_initial_condition(comp_v_array{5});

% Set the initial condition of the assembly
obj.coordinate          = obj.coordinate.set_initial_condition(q_0);
obj.velocity            = obj.velocity.set_initial_condition(v_0);
obj.velocity_reduced    = obj.velocity_reduced.set_initial_condition(v_hat_0);


end