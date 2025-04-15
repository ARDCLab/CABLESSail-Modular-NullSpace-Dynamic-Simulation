function [Force_sail, Torque_sail] = Compute_Sail_Force_SAIL_OBJECT(obj, p_ba, q1, q2, q3, q4)
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
[C_ba] = p2DCM(p_ba);
% THIS SHOULD BE C_AB, IS THERE A BUG SOMEWHERE ELSE
%C_ba = C_ab';
% THIS MIGHT BE WRONG???

% Transform sun into the body frame
% Assuming the be frame is aligne with the a grame at t=0
% Use a 2-rotation to set the SIA in the plane of booms 1 and 3
% SIA and clock describe the initial relative position of the sun
Sun_dir_sun_frame = [0;0;1];
SIA = 17; % degrees
clock_angle = 0;

C_a_sun = [cosd(SIA) 0 sind(SIA);
          0 1 0;
          -sind(SIA) 0 cosd(SIA)];
% assuming C_ab = eye(3) at t=0
Sun_dir_b = C_ba * C_a_sun*Sun_dir_sun_frame;

% Compute force and torque from SRP
[Force_sail,Torque_sail] = obj.sailMembraneTroque(Sun_dir_b, r_tcb_b(3), r_tdb_b(3), r_teb_b(3), r_tfb_b(3));

Force_sail = reshape(Force_sail, [3,1]);
Torque_sail = reshape(Torque_sail, [3,1]);
end












