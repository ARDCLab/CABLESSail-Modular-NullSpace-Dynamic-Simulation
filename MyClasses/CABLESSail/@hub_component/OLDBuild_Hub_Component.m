function comp = Build_Hub_Component(varagin)
% File to assemble Hub EoM constants
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: September 25 2023
%-------------------------------------------------------------------------%
if nargin==1
    comp_number = varagin;
else
    comp_number = 0;
end

% Material constants -----------------------------------------------------%
LoadConstants
% mass
hub.m_hub = 50; % kg
% height
hub.h = 1; % meter/Users/keeganbunker/CABLESail_sim/Figures_05-Dec-2023 09_32_24/1.png

% radius
hub.R = 1; % meter

m_sail = 50; % kg
a = 29.5*sqrt(2);
I_sail = m_sail*1/21*diag([a^2, a^2, 2*a^2]);
hub.m = m_sail + hub.m_hub;


% Moment of Inertia matrix
hub.I_Bbb_b = hub.m_hub*[1/12*hub.h^2 + 1/4*hub.R^2, 0, 0; ...
              0, 1/12*hub.h^2 + 1/4*hub.R^2, 0; ...
              0, 0, 1/2*hub.R^2]+I_sail;



% Mass and Stiffness Matrices --------------------------------------------%


% M mastrix for hub
hub.M = blkdiag(hub.m*eye(3), hub.I_Bbb_b);

% Define individual states
r_B_ba = state(3, 'Position',[0 0 0]');
v_B_ba = state(3, 'Velocity', [0 0 0]');
p_ba = state(9, 'Vectorized DCM', [1 0 0 0 1 0 0 0 1]');
w_ba = state(3, 'Angular Velocity', [0.0*const.deg2rad 0 0]');

% Assemble componet state vector and velocities
q_hub = state_vector({r_B_ba, p_ba});
vel_hub = state_vector({v_B_ba, w_ba});

% create component object
comp = component(comp_number, q_hub, vel_hub, @build_M_hub, @build_K_hub, @build_fnon_hub, @build_f_hub, @y2qv_and_q_dot_hub, hub, @build_Gamma_hub);

if nargin == 0
    tf = 60 ;
    dtN = tf*10;
    comp = comp.simulate_component(tf, dtN);
    %comp = comp.simulate_component(tf, 10);
    comp.state.plot;
    comp.compute_energies()
    comp.velocity.plot;
end


end % End function






