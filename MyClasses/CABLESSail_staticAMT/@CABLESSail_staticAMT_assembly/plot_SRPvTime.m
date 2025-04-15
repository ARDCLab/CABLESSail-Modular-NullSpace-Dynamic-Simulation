function plot_SRPvTime(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%
%
% ----------------------------------------------------------------------- %

% Bus cg location wrt inertial poit a
r_ba_a = obj.coordinate.states{1}.history;

% Boom Length
L_Bc = obj.components{2}.const.Boom_Length;
L_Bd = obj.components{3}.const.Boom_Length;
L_Be = obj.components{4}.const.Boom_Length;
L_Bf = obj.components{5}.const.Boom_Length;

% Boom tip deflections (point t) in the boom frames
r_tc_c = obj.coordinate.states{5}.deflection(obj.coordinate.states{5}.Length);
r_tc_c = [L_Bc*ones(1,length(r_tc_c)); r_tc_c];
r_td_d = obj.coordinate.states{8}.deflection(obj.coordinate.states{8}.Length);
r_td_d = [L_Bd*ones(1,length(r_td_d)); r_td_d];
r_te_e = obj.coordinate.states{11}.deflection(obj.coordinate.states{11}.Length);
r_te_e = [L_Be*ones(1,length(r_te_e)); r_te_e];
r_tf_f = obj.coordinate.states{14}.deflection(obj.coordinate.states{14}.Length);
r_tf_f = [L_Bf*ones(1,length(r_tf_f)); r_tf_f];

% Tip positions in Bus frame (sail body frame)
r_tc_b = obj.const.r_cb_b + obj.const.C_cb' * r_tc_c;
r_td_b = obj.const.r_db_b + obj.const.C_db' * r_td_d;
r_te_b = obj.const.r_eb_b + obj.const.C_eb' * r_te_e;
r_tf_b = obj.const.r_fb_b + obj.const.C_fb' * r_tf_f;

% Transform sun into the body frame -----------------------------------
% Assuming the be frame is aligne with the a grame at t=0
% Use a 2-rotation to set the SIA in the plane of booms 1 and 3
Sun_dir_sun_frame = [0;0;1];
SIA = 0; % degrees
% SIA IS WITH RESPECT TO THE NORMAL!!!!!!! NEED TO FIX
% CHANGE IN COMPONENT ASSEMPLY ODE
C_a_sun = [cosd(SIA) 0 sind(SIA);
          0 1 0;
          -sind(SIA) 0 cosd(SIA)];
% assuming C_ab = eye(3) at t=0
Sun_dir_a = C_a_sun*Sun_dir_sun_frame;

% Inialize force and torque matrices
Force_normal = zeros(size(r_ba_a));
Force_tangent = Force_normal;
Force_total = Force_normal;

Torque_normal = Force_normal;
Torque_tangent = Force_normal;
Torque_total = Force_normal;

% Loop through time histroy and compute torque and force
for lv11 = 1:length(r_ba_a)    
    % Form C_ab from time history
    p_ba_t = obj.coordinate.states{2}.history(:,lv11);
    [C_ba] = p2DCM(p_ba_t);
    C_ab = C_ba';
    % C_ab = eye(3);

    % Transform sun into the body frame
    Sun_dir_b = C_ab' * Sun_dir_a;
    %Sun_dir_b = Sun_dir_a;

    % Bus location
    Pos_val.hub = [0 0 0]';

    % Call rigid sail model force and torque computation function
    Pos_val.tip1 = r_tc_b(:,lv11);
    Pos_val.tip2 = r_td_b(:,lv11);
    [Force_t1, Torque_t1] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

    % Call rigid sail model force and torque computation function
    Pos_val.tip1 = r_td_b(:,lv11);
    Pos_val.tip2 = r_te_b(:,lv11);
    [Force_t2, Torque_t2] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

    % Call rigid sail model force and torque computation function
    Pos_val.tip1 = r_te_b(:,lv11);
    Pos_val.tip2 = r_tf_b(:,lv11);
    [Force_t3, Torque_t3] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);

    % Call rigid sail model force and torque computation function
    Pos_val.tip1 = r_tf_b(:,lv11);
    Pos_val.tip2 = r_tc_b(:,lv11);
    [Force_t4, Torque_t4] = RigidSailModel(Pos_val,Sun_dir_b,obj.const.Opt_prop);
    
    % for testing
    % [Torque_t1.total, Torque_t2.total , Torque_t3.total , Torque_t4.total ]
    % [r_tc_b(:,lv11), r_td_b(:,lv11), r_te_b(:,lv11), r_tf_b(:,lv11)]
    
    % Sum force and torque from all sail pieces
    Torque_total(:,lv11) = Torque_t1.total + Torque_t2.total + Torque_t3.total + Torque_t4.total;
    Force_total(:,lv11) = Force_t1.total + Force_t2.total + Force_t3.total + Force_t4.total;

    if ~isreal(Torque_total) || ~isreal(Force_total)
        disp("Complex Torque or Force")
    end
end

%% PLot Force results
% figure()
% hold on
% plot(obj.coordinate.states{1}.time, Torque_total(1,:) )
% plot(obj.coordinate.states{1}.time, Torque_total(2,:) )
% plot(obj.coordinate.states{1}.time, Torque_total(3,:) )
% legend('b1','b2','b3')
% ylabel("Generated Torque in Body Frame (N-m)")
% xlabel("Time (s)")
% grid on


%%

u1 = obj.coordinate.states{5}.deflection(obj.coordinate.states{5}.Length);
u2 = obj.coordinate.states{8}.deflection(obj.coordinate.states{5}.Length);
u3 = obj.coordinate.states{11}.deflection(obj.coordinate.states{5}.Length);
u4 = obj.coordinate.states{14}.deflection(obj.coordinate.states{5}.Length);

figure()
subplot(2,1,1)
plot( obj.coordinate.states{5}.time, u1(2,:) )
hold on
plot( obj.coordinate.states{5}.time, u2(2,:) )
plot( obj.coordinate.states{5}.time, u3(2,:) )
plot( obj.coordinate.states{5}.time, u4(2,:) )
xlim([obj.coordinate.states{5}.time(1), obj.coordinate.states{5}.time(end)]);
grid on
ylabel('Out of Plane Deflection $u_3$ (m)', 'interpreter','Latex');%, 'FontSize',8)
xlabel('Time (s)', 'interpreter','latex');%, 'FontSize',8)
leg = legend('$B_c$', '$B_d$', '$B_e$', '$B_f$', 'interpreter','Latex');%, 'FontSize',8);
leg.ItemTokenSize = [10,18];
%title(" Tip deflection of booms out of plane")

subplot(2,1,2)
hold on
plot(obj.coordinate.states{1}.time, Torque_total(1,:) )
plot(obj.coordinate.states{1}.time, Torque_total(2,:) )
plot(obj.coordinate.states{1}.time, Torque_total(3,:) )
xlim([obj.coordinate.states{5}.time(1), obj.coordinate.states{5}.time(end)]);
leg = legend('$b_1$','$b_2$','$b_3$');%, 'FontSize',8);
leg.ItemTokenSize = [10,18];
ylabel("Torque in Body Frame (N-m)");%, 'FontSize',8)
xlabel("Time (s)");%, 'FontSize',8)
grid on



end