function [Force_sail, Torque_sail] = Compute_Sail_Force_FlatePlate(obj, p_ba, q1, q2, q3, q4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% ----------------------------------------------------------------------- %

% Boom tip deflections (point t) in the boom frames
r_tc_c = [(obj.components{2}.const.Length); obj.components{2}.const.PSI_FH(obj.components{2}.const.Length)*q1];
r_td_d = [(obj.components{3}.const.Length); obj.components{3}.const.PSI_FH(obj.components{3}.const.Length)*q2];
r_te_e = [(obj.components{4}.const.Length); obj.components{4}.const.PSI_FH(obj.components{4}.const.Length)*q3];
r_tf_f = [(obj.components{5}.const.Length); obj.components{5}.const.PSI_FH(obj.components{5}.const.Length)*q4];

% Tip positions in Bus frame (sail body frame)
r_tcb_b = obj.const.r_cb_b + obj.const.C_cb' * r_tc_c;
r_tdb_b = obj.const.r_db_b + obj.const.C_db' * r_td_d;
r_teb_b = obj.const.r_eb_b + obj.const.C_eb' * r_te_e;
r_tfb_b = obj.const.r_fb_b + obj.const.C_fb' * r_tf_f;

% Sun direction in the inertial frame

% Bus location
Pos_val.hub = [0 0 0]';

% Form C_ab from time history
%[C_ab] = p2DCM(p_ba);
[C_ba] = p2DCM(p_ba);
C_ab = C_ba';

% Transform sun into the body frame
% Assuming the be frame is aligne with the a grame at t=0
% Use a 2-rotation to set the SIA in the plane of booms 1 and 3
Sun_dir_sun_frame = [0;0;1];
SIA = 17; % degrees
% CHANGE IN SRP PLOT SCRIPT
C_a_sun = [cosd(SIA) 0 sind(SIA);
          0 1 0;
          -sind(SIA) 0 cosd(SIA)];
% assuming C_ab = eye(3) at t=0
Sun_dir_b = C_ab' * C_a_sun*Sun_dir_sun_frame;

% Compute force and torque from SRP
% Sail 1
Pos_val.tip1 = r_tcb_b;
Pos_val.tip2 = r_tdb_b;
[Force_sail1,Torque_sail1] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

% Sail 2
Pos_val.tip1 = r_tdb_b;
Pos_val.tip2 = r_teb_b;
[Force_sail2,Torque_sail2] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

% Sail 3
Pos_val.tip1 = r_teb_b;
Pos_val.tip2 = r_tfb_b;
[Force_sail3,Torque_sail3] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

% Sail 4
Pos_val.tip1 = r_tfb_b;
Pos_val.tip2 = r_tcb_b;
[Force_sail4,Torque_sail4] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

% Total force and torques
Force_sail = Force_sail1.total + Force_sail2.total + Force_sail3.total + Force_sail4.total;
Torque_sail = Torque_sail1.total + Torque_sail2.total + Torque_sail3.total + Torque_sail4.total;

end