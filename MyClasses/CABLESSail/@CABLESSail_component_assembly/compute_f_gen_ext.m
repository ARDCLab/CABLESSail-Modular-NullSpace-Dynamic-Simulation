function f_gen_ext = compute_f_gen_ext(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% ----------------------------------------------------------------------- %

% Get the current elastic coordinates of each boom (Bc, Bd, Be, Bf)
p_ba     = comp_coords_array{1}{2};
q_Bc_eps = comp_coords_array{2}{3} ;
q_Bd_eps = comp_coords_array{3}{3} ;
q_Be_eps = comp_coords_array{4}{3} ;
q_Bf_eps = comp_coords_array{5}{3} ;

% Compute the generalized forces from the cables 
[f_gen_Bc] = obj.components{2}.compute_cable_generalized_forces(t, q_Bc_eps);
[f_gen_Bd] = obj.components{3}.compute_cable_generalized_forces(t, q_Bd_eps);
[f_gen_Be] = obj.components{4}.compute_cable_generalized_forces(t, q_Be_eps);
[f_gen_Bf] = obj.components{5}.compute_cable_generalized_forces(t, q_Bf_eps);

% Compute the generalized forces from the solar sail
%[F_gen_r, F_gen_w] = obj.compute_SS_generalized_forces_PLATE(p_ba,  q_Bc_eps, q_Bd_eps, q_Be_eps, q_Bf_eps);
[Force_sail, Torque_sail] = obj.Compute_Sail_Force_FlatePlate(p_ba, q_Bc_eps, q_Bd_eps, q_Be_eps, q_Bf_eps);
F_gen_r = Force_sail;
F_gen_w = Torque_sail;

% if t> 60*60
%     disp('breaking now');
% end

% Turn things off -----------------------
    F_gen_r = zeros(3,1);
    F_gen_w = zeros(3,1);
% ---------------------------------------

f_gen_ext = [F_gen_r; F_gen_w; f_gen_Bc; f_gen_Bd; f_gen_Be; f_gen_Bf];

end