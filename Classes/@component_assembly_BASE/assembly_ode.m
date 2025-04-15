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

% Extract state
% Decompose the state vector and compute Upsilon, Upsilon_dot
[coords, vels, vels_hat, q, v, v_hat, Upsilon, Upsilon_dot] = obj.decomepose_state(y);

% Create cell arrays of the component coordinates and velocites
[comp_coords_array, comp_vels_array, comp_q_array, comp_v_array] = obj.decompoese_component_states(coords, vels, q, v);
    

% Current state matrices
K           = obj.compute_K(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
M           = obj.compute_M(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
f_non       = obj.compute_f_non(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
f_gen       = obj.compute_f_gen(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
gamma       = obj.compute_Gamma(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
f_gen_ext   = obj.compute_f_gen_ext(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
%D           = obj.compute_D(t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);

% ----------------------------------------
% Checking the actual damping of the system
% ----------------------------------------
D = blkdiag(eye(3), eye(3), 3*10e-6*eye(6*4));
%D = zeros(30,30);

% % compute nat-freqs and damping ration from the linear system
% a11 = zeros(96,96);
% a12 = gamma*Upsilon;
% a21 = -inv(Upsilon'* M * Upsilon)*Upsilon'*gamma'*K;
% a22 = -inv(Upsilon'* M * Upsilon)*D;
% A_linear = [a11, a12;
%             a21, zeros(30,30)];
% A_linear_damp = [a11, a12;
%             a21, a22];
% 
% % Compute the natural frequenices in hz
% w_n2 = abs(eig(A_linear));
% w_nHZ = sqrt(w_n2)/(2*pi);
% 
% % Compute the damping ratio
% zeta = -cos(angle(eig(A_linear_damp)));
% sort(zeta(zeta>0))

% ----------------------------------------
% ----------------------------------------

% if t> 60*60
%     disp('breaking now');
% end

% Update state
v_dot = inv(Upsilon'* M * Upsilon) * ( Upsilon'*(gamma'*f_gen - M*Upsilon_dot*v_hat - gamma'*K*q - f_non) - D*v_hat + f_gen_ext );
dydt = [gamma*Upsilon*v_hat; v_dot];


end



























