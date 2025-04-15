function [Force,Torque] = RigidSailModel(Pos_val,Sun_dir,Opt_prop)
%RigidSailModel - Computes SRP forces and torques due to rigid triangular sail
%
% Syntax: [Force,Torque] = RigidSailModel(Pos_val,Sun_dir,Opt_prop)
%
% Inputs:
%   Pos_val  - Structure that includes position information of the corners 
%              of the sail quadrant in the body frame (assumed that
%              positive 3 axis of body frame is nominally pointed at Sun
%              when SIA = 0 deg)
%              Pos_val.hub  - position of the hub attachment point (m)
%              Pos_val.tip1 - position of the first boom-tip attachment
%                             point (m)
%              Pos_val.tip2 - position of the second boom-tip attachment
%                             point (m)
%   Sun_dir  - Vector pointing in the direction of the Sun resolved in the
%              body frame
%   Opt_prop - Structure containing the optical properties of the sail
%              (values in parentheses are for Solar Cruiser)
%              Opt_prop.P       - solar radiation (4.53914e-6 N/m^2 at 1AU)
%              Opt_prop.r_tilde - reflection coefficient (0.91)
%              Opt_prop.s       - fraction of specular reflection
%                                 coefficient (0.94)
%              Opt_prop.Bf      - front non-Lambertian coefficient (0.79)
%              Opt_prop.Bb      - back non-Lambertian coefficient (0.67)
%              Opt_prop.ef      - front emissivity (0.025)
%              Opt_prop.eb      - back emissivity (0.27)
% Outputs:
%   Force    - Structure with SRP forces
%              Force.total   - total force from normal and tangent
%              Force.normal  - normal force only
%              Force.tangent - tangent force only
%   Torque   - Structure with SRP torques
%              Torque.total   - total torque from normal and tangent forces
%              Torque.normal  - torque from only normal SRP force
%              Torque.tangent - torque from only tangent SRP force
%
% Author: Ryan Caverly
% November 2023; Last revision: 13-Nov-2023

%------------- BEGIN CODE --------------

% Reformat data into column vectors
% Position of the hub resolved in the body frame
Pos_val.hub = [Pos_val.hub(1);Pos_val.hub(2);Pos_val.hub(3)];

% Position of the tips resolved in the body frame
Pos_val.tip1 = [Pos_val.tip1(1);Pos_val.tip1(2);Pos_val.tip1(3)];
Pos_val.tip2 = [Pos_val.tip2(1);Pos_val.tip2(2);Pos_val.tip2(3)];
Sun_dir = [Sun_dir(1);Sun_dir(2);Sun_dir(3)]/norm(Sun_dir);

% Compute position vectors for each side of the sail
vec1 = Pos_val.tip1 -  Pos_val.hub; % Position vector of 1st sail side
vec2 = Pos_val.tip2 -  Pos_val.hub; % Position vector of 2nd sail side
sail_normal = cross(vec1,vec2)/norm(cross(vec1,vec2)); % Normal vector to sail
if sail_normal(3) < 0 % Ensure that normal vector is pointing in positive direction
    sail_normal = -sail_normal;
end

% Compute sail area and the center of the sail
sail_area = norm(cross(vec1,vec2))/2; % Area of sail
sail_center =  (Pos_val.hub + Pos_val.tip1 + Pos_val.tip2)/3; % Center of sail

% Vector tangent to sail
if norm(Sun_dir-sail_normal) == 0
    sail_tangent = sail_normal;
else
    sail_tangent = (Sun_dir-sail_normal)/norm(Sun_dir-sail_normal); 
end

% Compute sun incidence angle of the sail
Local_SIA = acos(sail_normal'*Sun_dir);
% Above sometimes returns a complex value
if ~isreal(Local_SIA)
    disp("Complex SIA calculated: ")
    disp(Local_SIA)
    Local_SIA = real(Local_SIA);
end
cos_theta = cos(Local_SIA);
sin_theta = sin(Local_SIA);

% magnitude of SRP in normal direction
force_mag_normal = Opt_prop.P*sail_area*(1+Opt_prop.r_tilde*Opt_prop.s)*cos_theta^2 ...
    + Opt_prop.P*sail_area*Opt_prop.Bf*(1-Opt_prop.s)*Opt_prop.r_tilde*cos_theta ...
    + Opt_prop.P*sail_area*(1-Opt_prop.r_tilde)*(Opt_prop.ef*Opt_prop.Bf ...
    -Opt_prop.eb*Opt_prop.Bb)/(Opt_prop.ef+Opt_prop.eb)*cos_theta; 

% magnitude of SRP in tangent direction
force_mag_tangent = Opt_prop.P*sail_area*(1-Opt_prop.r_tilde*Opt_prop.s)*cos_theta*sin_theta; 

% Compute forces
Force.normal = -force_mag_normal*sail_normal;
Force.tangent = - force_mag_tangent*sail_tangent;
Force.total = Force.normal + Force.tangent;

% Compute torques
Torque.normal = cross(sail_center,Force.normal);
Torque.tangent = cross(sail_center,Force.tangent);
Torque.total = Torque.normal + Torque.tangent;
