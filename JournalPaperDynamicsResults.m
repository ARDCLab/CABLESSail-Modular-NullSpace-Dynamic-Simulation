% Keegan R Bunker
% PhD Student
% University of Minnesota
% PI: Dr. Ryan Caverly
% Initial date: 8/1/2024
% Last edits made:
% 
% This is the main script to run the tests cases for the AIAA Journal Paper
% extension of the scitech paper.
%
% Two test cases are ran in this script
%   1. A CABLESSail'ed solar cruiser spacecraft with one boom being slowly
%   actuated to a constant tip deflection
%
%   2. A solar-Cruiser-like spacecraft with a fully actuated AMT shifting
%   the center of mass of the bus-AMT. 
%
% Both of these test cases are subjected to the same conditions.
%   - Zero initial velocity and acceleration of all states
%   - Flat plate membrane shape for each quadrant
%   - Non zero sun incidence angle (SIA = 17, clock = 0? check other work)
%   - Over the course of hours (?) examine how the attitude changes over
%   time due to the SRP torques
%

%% Workspace and preliminaries -------------------------------------------

% Add the classes and utilities folder of the simulation to the workspace
addpath(genpath([pwd, '/Utilities']))
addpath(genpath([pwd, '/Classes']))
set(0,'DefaultFigureWindowStyle','docked')

clear all

%% Create a non flat sail object for the tests
addpath(genpath([pwd, '/MyClasses/Sail_Shape_Object']))
% Sail Parameters
L_boom                      = 30;% Boom Length
n_nodes                     = 30;% Number of nodes
n_sailShape                 = 2;
basisFunction               = @(x,y) BasisGauvainTyler(x, y, L_boom);
CABLESSail_membrane_billow = SailMesh(...
    L_boom, ...
    n_nodes, ...
    {basisFunction, basisFunction, basisFunction, basisFunction}, ...
    {-0.05, -0.05, -0.05, -0.05 });

CABLESSail_membrane_flat = SailMesh(...
    L_boom, ...
    2, ...
    {@(x,y) 0, @(x,y) 0, @(x,y) 0, @(x,y) 0}, ...
    {1, 1, 1, 1 });


%% Test Case 1: CABLESSail MEMBRANE ---------------------------------------
% ----------------------------------------------------------------------- %

% Add the folder containing the classes I created to the workspace
addpath(genpath([pwd, '/MyClasses/CABLESSail_staticAMT']))

% Create the components for this assembly
Bb_hub = bus_component();
Bc_boom1 = boom_component();
Bd_boom2 = boom_component();
Be_boom3 = boom_component();
Bf_boom4 = boom_component();
Bg_AMT = AMT_component();
components = {Bb_hub, Bc_boom1, Bd_boom2, Be_boom3, Bf_boom4, Bg_AMT};

% Create the vector of reduced states
v_ba_a = state(3, 'Velocity');
w_ba_b = state(3, 'Angular Velocity');
q1_dot = state(6, 'd/dt elastic_coordinates');
q2_dot = state(6, 'd/dt elastic_coordinates');
q3_dot = state(6, 'd/dt elastic_coordinates');
q4_dot = state(6, 'd/dt elastic_coordinates');
v_hat = state_vector({v_ba_a, w_ba_b, q1_dot, q2_dot, q3_dot, q4_dot});

% Form the component assembly
assembly = CABLESSail_staticAMT_assembly(components, v_hat);

% Set the SRP torque function to the flat mesh
assembly.sailMembraneTroque = @(Sun_dir_b, eps1, eps2, eps3, eps4) CABLESSail_membrane_flat.computeSRPForceandTorque(Sun_dir_b, eps1, eps2, eps3, eps4);

% Set the initial condition of generalized coordinates
r_ba_a0 = [0 0 0]';
p_ba0 = p2DCM(eye(3));
q10 = [0 0 0 0 0 0]';
q20 = [0 0 0 0 0 0]';
q30 = [0 0 0 0 0 0]';
q40 = [0 0 0 0 0 0]';

% Set the initial conditions of the generalized velocities
deg2rad = pi/180;
v_ba_a0 = [0 0 0]';
w_ba_b0 = [0 0 0]';
%w_ba_b0 = [0 1*deg2rad 1*deg2rad]';
q1_dot0 = [0 0 0 0 0 0]';
q2_dot0 = [0 0 0 0 0 0]';
q3_dot0 = [0 0 0 0 0 0]';
q4_dot0 = [0 0 0 0 0 0]';

% Assemble initial state variables
coords_reduced_0    = {r_ba_a0, p_ba0, q10, q20, q30, q40};
q_hat_0             = [r_ba_a0; p_ba0; q10; q20; q30; q40];
vels_reduced_0      = {v_ba_a0, w_ba_b0, q1_dot0, q2_dot0, q3_dot0, q4_dot0};
v_hat_0             = [v_ba_a0; w_ba_b0; q1_dot0; q2_dot0; q3_dot0; q4_dot0];
assembly            = assembly.set_initial_conditions(coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0);
    % Update the set initial conditions method so you don't have to pass
    % the same data as an array and a cell array

% Set the cable tension
%assembly.components{2}.T3 = @(t) 0.1; % Newtons
%maxTension = 0.04; %Newtons
maxTension = 4.00; %Newtons
maxTimeSecs = 60; % Seconds. One hour ramp time
assembly.components{2}.T3 = @(t) TensionRamp(t, maxTension, maxTimeSecs); % Newtons

% Simulate test 1
tf = 30*60; %*2
%tf = 60;
dtN = tf/10;
assembly = assembly.simulate(tf, dtN);

%% Plot cablessail test case
% assembly.coordinate.states{1}.plot;
assembly.coordinate.states{2}.plot;
assembly.coordinate.states{5}.plot;
% assembly.coordinate.states{8}.plot;
% assembly.coordinate.states{11}.plot;
% assembly.coordinate.states{14}.plot;
assembly.plot_SRPvTime();

%% Test Case 2: CABLESSail FLAT PLATE -------------------------------------
% ----------------------------------------------------------------------- %
% clear assembly
% % Add the folder containing the classes I created to the workspace
% addpath(genpath([pwd, '/MyClasses/CABLESSail_staticAMT']))
% 
% % Create the components for this assembly
% Bb_hub = bus_component();
% Bc_boom1 = boom_component();
% Bd_boom2 = boom_component();
% Be_boom3 = boom_component();
% Bf_boom4 = boom_component();
% Bg_AMT = AMT_component();
% components = {Bb_hub, Bc_boom1, Bd_boom2, Be_boom3, Bf_boom4, Bg_AMT};
% 
% % Create the vector of reduced states
% v_ba_a = state(3, 'Velocity');
% w_ba_b = state(3, 'Angular Velocity');
% q1_dot = state(6, 'd/dt elastic_coordinates');
% q2_dot = state(6, 'd/dt elastic_coordinates');
% q3_dot = state(6, 'd/dt elastic_coordinates');
% q4_dot = state(6, 'd/dt elastic_coordinates');
% v_hat = state_vector({v_ba_a, w_ba_b, q1_dot, q2_dot, q3_dot, q4_dot});
% 
% % Form the component assembly
% assembly = CABLESSail_staticAMT_assembly(components, v_hat);
% 
% % Set the SRP torque function to the flat mesh
% assembly.sailMembraneTroque = @(Sun_dir_b, eps1, eps2, eps3, eps4) CABLESSail_membrane_billow.computeSRPForceandTorque(Sun_dir_b, eps1, eps2, eps3, eps4);
% 
% % Set the initial condition of generalized coordinates
% r_ba_a0 = [0 0 0]';
% p_ba0 = p2DCM(eye(3));
% q10 = [0 0 0 0 0 0]';
% q20 = [0 0 0 0 0 0]';
% q30 = [0 0 0 0 0 0]';
% q40 = [0 0 0 0 0 0]';
% 
% % Set the initial conditions of the generalized velocities
% deg2rad = pi/180;
% v_ba_a0 = [0 0 0]';
% w_ba_b0 = [0 0 0]';
% %w_ba_b0 = [0 1*deg2rad 1*deg2rad]';
% q1_dot0 = [0 0 0 0 0 0]';
% q2_dot0 = [0 0 0 0 0 0]';
% q3_dot0 = [0 0 0 0 0 0]';
% q4_dot0 = [0 0 0 0 0 0]';
% 
% % Assemble initial state variables
% coords_reduced_0    = {r_ba_a0, p_ba0, q10, q20, q30, q40};
% q_hat_0             = [r_ba_a0; p_ba0; q10; q20; q30; q40];
% vels_reduced_0      = {v_ba_a0, w_ba_b0, q1_dot0, q2_dot0, q3_dot0, q4_dot0};
% v_hat_0             = [v_ba_a0; w_ba_b0; q1_dot0; q2_dot0; q3_dot0; q4_dot0];
% assembly            = assembly.set_initial_conditions(coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0);
%     % Update the set initial conditions method so you don't have to pass
%     % the same data as an array and a cell array
% 
% % Set the cable tension
% %assembly.components{2}.T3 = @(t) 0.1; % Newtons
% maxTension = 0.04; %Newtons
% maxTimeSecs = 2*60*60; % Seconds. One hour ramp time
% assembly.components{2}.T3 = @(t) TensionRamp(t, maxTension, maxTimeSecs); % Newtons
% 
% % Simulate test 1
% tf = 60*60;
% dtN = tf/10;
% assembly = assembly.simulate(tf, dtN);

%% Plot cablessail test case
% % assembly.coordinate.states{1}.plot;
% assembly.coordinate.states{2}.plot;
% assembly.coordinate.states{5}.plot;
% % assembly.coordinate.states{8}.plot;
% % assembly.coordinate.states{11}.plot;
% % assembly.coordinate.states{14}.plot;
% assembly.plot_SRPvTime();

%% Test Case Two: AMT
% ----------------------------------------------------------------------- %
clear
clc
% Add the folder containing the classes I created to the workspace
addpath(genpath([pwd, '/MyClasses/CABLESSail_staticAMT']))

% Create the components for this assembly
Bb_hub = bus_component();
Bc_boom1 = boom_component();
Bd_boom2 = boom_component();
Be_boom3 = boom_component();
Bf_boom4 = boom_component();
Bg_AMT = AMT_component();
components = {Bb_hub, Bc_boom1, Bd_boom2, Be_boom3, Bf_boom4, Bg_AMT};

% Create the vector of reduced states
v_ba_a = state(3, 'Velocity');
w_ba_b = state(3, 'Angular Velocity');
q1_dot = state(6, 'd/dt elastic_coordinates');
q2_dot = state(6, 'd/dt elastic_coordinates');
q3_dot = state(6, 'd/dt elastic_coordinates');
q4_dot = state(6, 'd/dt elastic_coordinates');
v_hat = state_vector({v_ba_a, w_ba_b, q1_dot, q2_dot, q3_dot, q4_dot});

% Form the component assembly
assembly = CABLESSail_staticAMT_assembly(components, v_hat);

% Set the initial condition of generalized coordinates
r_ba_a0 = [0 0 0]';
p_ba0 = p2DCM(eye(3));
q10 = [0 0 0 0 0 0]';
q20 = [0 0 0 0 0 0]';
q30 = [0 0 0 0 0 0]';
q40 = [0 0 0 0 0 0]';

% Set the initial conditions of the generalized velocities
deg2rad = pi/180;
v_ba_a0 = [0 0 0]';
w_ba_b0 = [0 0 0]';
%w_ba_b0 = [0 1*deg2rad 1*deg2rad]';
q1_dot0 = [0 0 0 0 0 0]';
q2_dot0 = [0 0 0 0 0 0]';
q3_dot0 = [0 0 0 0 0 0]';
q4_dot0 = [0 0 0 0 0 0]';

% Assemble initial state variables
coords_reduced_0    = {r_ba_a0, p_ba0, q10, q20, q30, q40};
q_hat_0             = [r_ba_a0; p_ba0; q10; q20; q30; q40];
vels_reduced_0      = {v_ba_a0, w_ba_b0, q1_dot0, q2_dot0, q3_dot0, q4_dot0};
v_hat_0             = [v_ba_a0; w_ba_b0; q1_dot0; q2_dot0; q3_dot0; q4_dot0];
assembly            = assembly.set_initial_conditions(coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0);
    % Update the set initial conditions method so you don't have to pass
    % the same data as an array and a cell array

% Change the AMT position and set tension to 0
assembly.const.r_gb_b = [0.30; 0; assembly.const.r_gb_b(3)];
assembly.components{2}.T3 = @(t) 0.0; % Newtons

% Simulate
tf = 30*60;
dtN = tf/10;
assembly = assembly.simulate(tf, dtN);

%% Plot AMT case
% assembly.coordinate.states{1}.plot;
assembly.coordinate.states{2}.plot;
% assembly.coordinate.states{5}.plot;
% assembly.coordinate.states{8}.plot;
% assembly.coordinate.states{11}.plot;
% assembly.coordinate.states{14}.plot;

%% local functions 
% ----------------------------------------------------------------------- %

function T = TensionRamp(t, maxT, maxTimeSecs)
    if t<=maxTimeSecs
        T =t/maxTimeSecs * maxT;
    else
        T = maxT;
    end
end












