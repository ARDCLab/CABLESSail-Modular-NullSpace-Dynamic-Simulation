function [Force,Torque] = MeshSailModel(X, Y, Z, T, Sun_dir,Opt_prop)
% MeshSailModel - Computes the SRP forces and torques due to a deformed
% sail made up of a triangular mesh.
%
% Inputs:
%   X, Y, Z   - Column vectors of x, y, and z coordinates respectfully of
%               all mesh vertex points. Index i of each vector creates the
%               the vector defining the position of the i'th index point
%               such that r_i = [X(i) Y(i) Z(i)]'.
%
%   T         - The nx3 2D delaunay triangulation grid corrsponding to
%               the triangular mesh defined by X, Y, and Z.
%
%   Sun_dir   - The unit vector pointing in the direction of the sun
%               resolved in the same frame the mesh (X, Y, Z).
%
%   Opt_prop  - Structure containing the optical properties of the sail
%               (values in parentheses are for Solar Cruiser)
%               Opt_prop.P       - solar radiation (4.53914e-6 N/m^2 at 1AU)
%               Opt_prop.r_tilde - reflection coefficient (0.91)
%               Opt_prop.s       - fraction of specular reflection
%                                 coefficient (0.94)
%               Opt_prop.Bf      - front non-Lambertian coefficient (0.79)
%               Opt_prop.Bb      - back non-Lambertian coefficient (0.67)
%               Opt_prop.ef      - front emissivity (0.025)
%               Opt_prop.eb      - back emissivity (0.27)
%
% Outputs:
%   Force     - Structure with SRP forces
% 
%               Force.total   - total force from normal and tangent
%               Force.normal  - normal force only
%               Force.tangent - tangent force only
%   Torque    - Structure with SRP torques
%               Torque.total   - total torque from normal and tangent forces
%               Torque.normal  - torque from only normal SRP force
%               Torque.tangent - torque from only tangent SRP force
%
% Author: Keegan Bunker, Ryan Caverly
% March 2024; Last revision:
% ----------------------------------------------------------------------- %

% Preallocate matrices for speed
face_force_normal   = zeros(length(X),3);
face_force_tangent  = zeros(length(X),3);
face_torque_normal  = zeros(length(X),3);
face_torque_tangent = zeros(length(X),3);

for lv1=1:length(T) % Loop through all faces within mesh

    % Extract the vertices of the current triangle based on the T
    point1 = [X(T(lv1, 1)) Y(T(lv1, 1)) Z(T(lv1, 1))]';
    point2 = [X(T(lv1, 2)) Y(T(lv1, 2)) Z(T(lv1, 2))]';
    point3 = [X(T(lv1, 3)) Y(T(lv1, 3)) Z(T(lv1, 3))]';
                
    % Create unit vectors to define the surface elemnt
    vec1 = point2 - point1; % Position vector of 1st face side
    vec2 = point3 - point1;
    vec3 = cross(vec1,vec2)/norm(cross(vec1,vec2)); % Normal unit vector to face
    if vec3(3) < 0 % Ensure that normal vector is pointing in positive direction
       vec3 = -vec3;
    end
    
    %face_normal(lv1,:) = vec3; % Save normal vector
    area = norm(cross(vec1,vec2))/2; % Area of face
    face_center = (point1 + point2 + point3)/3; % Position of center of face
    
    if norm(Sun_dir-vec3) == 0
       face_tangent = vec3;
    else
       face_tangent = (Sun_dir-vec3)/norm(Sun_dir-vec3); % vector tangent to face
    end
    
    % Compute sun incidence angle of the sail
    face_angle = acos(vec3'*Sun_dir); % angle between sun vector and face normal
    cos_face_angle = cos(face_angle);

    % magnitude of SRP in normal direction
    force_mag_normal = Opt_prop.P*area*(1+Opt_prop.r_tilde*Opt_prop.s)*cos_face_angle^2 ...
    + Opt_prop.P*area*Opt_prop.Bf*(1-Opt_prop.s)*Opt_prop.r_tilde*cos_face_angle ...
    + Opt_prop.P*area*(1-Opt_prop.r_tilde)*(Opt_prop.ef*Opt_prop.Bf ...
    -Opt_prop.eb*Opt_prop.Bb)/(Opt_prop.ef+Opt_prop.eb)*cos_face_angle; 

    % magnitude of SRP in tangent direction
    force_mag_tangent = Opt_prop.P*area*(1-Opt_prop.r_tilde*Opt_prop.s)*cos_face_angle*sin(face_angle); % magnitude of SRP in tangent direction

    % Create force vectors with magnitude and unit vector of face
    face_force_normal(lv1,:) = -force_mag_normal*vec3';
    face_force_tangent(lv1,:) = - force_mag_tangent*face_tangent'; % total SRP force on the face

    face_torque_normal(lv1,:) = cross(face_center,face_force_normal(lv1,:)); % total SRP torque on the face
    face_torque_tangent(lv1,:) = cross(face_center,face_force_tangent(lv1,:)); % total SRP torque on the face
end

Force = sum(face_force_normal) + sum(face_force_tangent); % total SRP force summed across all faces
Torque = sum(face_torque_normal) + sum(face_torque_tangent);% total SRP torque summed across all faces


end