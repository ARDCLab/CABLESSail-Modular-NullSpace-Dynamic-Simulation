function [coords, q] = q_hat2q(obj, coords_reduced, q_hat)
% q_hat2q Computes the expanded set of coordinates from the reduced set of
% coordinates.
%   This enforces the holonomic constraints at the coordinate level. Since
%   they are holonomic and enforced on the velocities, this constraint only
%   needs to be enforced at the initial condition. If not enforced, the
%   coordinates are simply off by a constant offset.
%
% Inputs
%   obj             - Parent class object
%   coords_reduced  - Cell array of seperate vector states that make up the
%                     full state vector q_hat 
%   q_hat           - Coordinate state in column
%   vector form
% 
% Outputs
%   coords - cell array of seperate state vectors that make up q
%   q      - Coordinate state in column vector form
%
% ----------------------------------------------------------------------- %

r_ba_a = coords_reduced{1};
p_ba = coords_reduced{2};
q1 = coords_reduced{3};
q2 = coords_reduced{4};
q3 = coords_reduced{5};
q4 = coords_reduced{6};

% Get DCM
C_ba = p2DCM(p_ba);
C_ab = C_ba';

% Boom components ------------------------------
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

% Resolve boom locations in the inerital frame
r_cb_a = C_ab* r_cb_b;
r_db_a = C_ab* r_db_b;
r_eb_a = C_ab* r_eb_b;
r_fb_a = C_ab* r_fb_b;

r_ca_a = r_cb_a  + r_ba_a;
r_da_a = r_db_a  + r_ba_a;
r_ea_a = r_eb_a  + r_ba_a;
r_fa_a = r_fb_a  + r_ba_a;

% Compute DCMs of booms
C_ca = C_cb*C_ba;
C_da = C_db*C_ba;
C_ea = C_eb*C_ba;
C_fa = C_fb*C_ba;

% Decompose boom DCM's
[p_ca] = DCM2p(C_ca);
[p_da] = DCM2p(C_da);
[p_ea] = DCM2p(C_ea);
[p_fa] = DCM2p(C_fa);

% AMT component ----------------------------
% Get the AMT location
r_gb_b = obj.const.r_gb_b;
r_gb_a = C_ab*r_gb_b;
r_ga_a = r_gb_a + r_ba_a;
C_ga = C_ba;
[p_ga] = DCM2p(C_ga);

% Form state vectors -----------------------
q_Bb_0 = [r_ba_a; p_ba];
q_Bc_0 = [r_ca_a; p_ca; q1];
q_Bd_0 = [r_da_a; p_da; q2];
q_Be_0 = [r_ea_a; p_ea; q3];
q_Bf_0 = [r_fa_a; p_fa; q4];
q_Bg_0 = [r_ga_a; p_ga];

% Compute all states and velocities
q = [q_Bb_0; q_Bc_0; q_Bd_0; q_Be_0; q_Bf_0; q_Bg_0];
coords = {r_ba_a, p_ba, ...
        r_ca_a, p_ca, q1, ...
        r_da_a, p_da, q2, ...
        r_ea_a, p_ea, q3, ...
        r_fa_a, p_fa, q4, ...
        r_ga_a, p_ga};

end