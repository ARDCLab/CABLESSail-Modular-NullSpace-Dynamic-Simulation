% Main cable sail simulator
% Add folders to workspace ---------------------------------------------- %
addpath(genpath([pwd, '/Utilities']))
addpath(genpath([pwd, '/Classes']))
set(0,'DefaultFigureWindowStyle','docked')
ForceLatexInterpreter


% Create all components
Bb_hub = hub_component();
Bc_boom1 = boom_component();
Bd_boom2 = boom_component();
Be_boom3 = boom_component();
Bf_boom4 = boom_component();
components = {Bb_hub, Bc_boom1, Bd_boom2, Be_boom3, Bf_boom4};

% Create the vector of reduced states
v_ba_a = state(3, 'Velocity');
w_ba_b = state(3, 'Angular Velocity');
q1_dot = state(6, 'd/dt elastic_coordinates');
q2_dot = state(6, 'd/dt elastic_coordinates');
q3_dot = state(6, 'd/dt elastic_coordinates');
q4_dot = state(6, 'd/dt elastic_coordinates');
v_hat = state_vector({v_ba_a, w_ba_b, q1_dot, q2_dot, q3_dot, q4_dot});

% Form the component assembly
assembly = CABLESSail_component_assembly(components, v_hat);

%% Simulate

% Set the initial conditions
r_ba_a0 = [0 0 0]';
v_ba_a0 = [0 0 0]';
p_ba0 = p2DCM(eye(3));
w_ba_b0 = [0 0 0]';
q10 = [0 0 0 0 0 0]';
q20 = [0 0 0 0 0 0]';
q30 = [0 0 0 0 0 0]';
q40 = [0 0 0 0 0 0]';
q1_dot0 = [0 0 0 0 0 0]';
q2_dot0 = [0 0 0 0 0 0]';
q3_dot0 = [0 0 0 0 0 0]';
q4_dot0 = [0 0 0 0 0 0]';

% Assemble state variables
coords_reduced_0    = {r_ba_a0, p_ba0, q10, q20, q30, q40};
q_hat_0             = [r_ba_a0; p_ba0; q10; q20; q30; q40];
vels_reduced_0      = {v_ba_a0, w_ba_b0, q1_dot0, q2_dot0, q3_dot0, q4_dot0};
v_hat_0             = [v_ba_a0; w_ba_b0; q1_dot0; q2_dot0; q3_dot0; q4_dot0];
assembly            = assembly.set_initial_conditions(coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0);
    % Update the set initial conditions method so you don't have to pass
    % the same data as an array and a cell array

% Set the cable tension
assembly.components{2}.T3 = 0.1; % Newtons

% Simulate
tf = 60*60;
dtN = tf/10;
assembly = assembly.simulate(tf, dtN);

%% Plot states

% % Plot positions
% assembly.coordinate.states{1}.plot;
% 
% % Plot attitude
% assembly.coordinate.states{2}.plot;

% Plot the first set of elastic coordinates
assembly.coordinate.states{5}.plot;
assembly.coordinate.states{8}.plot;
assembly.coordinate.states{11}.plot;
assembly.coordinate.states{14}.plot;

%% Confirm that expanded set of coordinates are changing as expect with the attitude evolving through time
% assembly.coordinate.states{1}.plot;
% assembly.coordinate.states{3}.plot;
% assembly.coordinate.states{6}.plot;
% assembly.coordinate.states{9}.plot;
% assembly.coordinate.states{12}.plot;


%% Boom Deflection and SRP torques
assembly.plot_SRPvTime()






