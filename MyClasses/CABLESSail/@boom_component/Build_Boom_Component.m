function comp = Build_Boom_Component(varagin)

% File to assemble boom EoM constants
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: August 14, 2023
%
%-------------------------------------------------------------------------%

if nargin==1
    comp_number = varagin;
else
    comp_number = 0;
end


% Material constants to match
% E_solar_cruiser = 80*10^9; % Pascals (Gpa = 1e9 = 10^9)
% I_2ma = 6.2e-8;% m^4


% Constants ------------------------------------------------------------- %
LoadConstants
boom.Length = 29.5; % meter
% boom.Modulus_Elasticity = 228*1e9; % Pascals (Gpa = 1e9 = 10^9)
boom.Modulus_Elasticity = 228*1e3; % Pascals (Gpa = 1e9 = 10^9)
%boom.Modulus_Elasticity = 228*1e5;
boom.Radius = 0.1;
boom.Plate_h = 0.1;
%boom.CrossSection_MOI = 4.99*10^(-10); % m^4
boom.CrossSection_MOI =pi/64*(2*boom.Radius)^4; % 2nd moment of area
boom.m = 3; % kg
boom.rho = boom.m/(pi*boom.Radius^2*boom.Length); 
boom.rho_L = boom.m/(boom.Length); 


% Basis Functions etc --------------------------------------------------- %
% Define the basis functions
syms x;
%psi_i = [(x./boom.Length).^2 (x./boom.Length).^3 (x./boom.Length).^4];% Pick some basis functions
psi_i = [(x/boom.Length).^2 (x/boom.Length).^3 (x/boom.Length).^4];% Pick some basis functions
% psi_i = [(x/boom.Length^(1/2)).^2 (x/boom.Length^(1/2)).^3 (x/boom.Length^(1/2)).^4];% Pick some basis functions
% psi_i = [(x/boom.Length^(2)).^2 (x/boom.Length^(2)).^3 (x/boom.Length^(2)).^4];% Pick some basis functions
n_psi = length(psi_i);
Psi = [psi_i,           zeros(1,n_psi); 
       zeros(1,n_psi),  psi_i ];
boom.PSI = Psi;
boom.PSI_FH = matlabFunction(Psi);
dpsi_dx2 = diff(Psi,x,2);
boom.psi_i_FH = matlabFunction(psi_i);
boom.dpsi_i_dx_FH = matlabFunction(diff(psi_i,x,1));

% Compute int Psi^T * Psi dx
P_B_i = double(int(boom.rho_L*(Psi'*Psi), x, 0, boom.Length));

% Spreader plate locaitons and psi matrix
boom.x_plates = linspace(0,29.5, 20);
boom.psi_X = boom.psi_i_FH(boom.x_plates');

% Positive def check ---------------------------------------------------- %
% Pdefint = [eye(3), [zeros(1,6); Psi];
%            [zeros(1,6); Psi]', [zeros(1,6); Psi]'*[zeros(1,6); Psi]];
% m_check = double(int(Pdefint, x, 0, boom.Length));

% Stiffness matrix ------------------------------------------------------ %
K_integrand = blkdiag(zeros(3), zeros(9), boom.Modulus_Elasticity*boom.CrossSection_MOI*(dpsi_dx2'*dpsi_dx2));
K = double(int(K_integrand, x, 0, boom.Length));
boom.K = K;
clear K_integrand;

% Mass Matrix ----------------------------------------------------------- %

% Symblic vector function used to check mass matrix
boom.r_dmi_i = [x; 0; 0];

% ----------------------------------------------------------------------- %

% Integrals for mass matrix
syms x y z;
Gint0 = crossr([x,y,z]')'*[zeros(1,6); Psi];
Gint1 = int(Gint0, z, -sqrt(boom.Radius^2 - y^2), sqrt(boom.Radius^2 - y^2));
Gint2 = int(Gint1, y, -boom.Radius, boom.Radius);
G_B_i = double(boom.rho*int(Gint2, x, 0, boom.Length));

% Construct H for booms
H_integrand = [zeros(1,2*n_psi); Psi];
H_B_i = double(int(boom.rho_L*H_integrand, x, [0, boom.Length]));
%PsiTransPsi2 = double(int(H_integrand*H_integrand', x, [0, boom.Length]));
clear H_integrand;

% Moment of inertia tensor of boom about the base
I_b_i = boom.m*[1/2*boom.Radius^2, 0, 0;
        0, 1/3*boom.Length^2 + 1/4*boom.Radius^2, 0;
        0, 0, 1/3*boom.Length^2 + 1/4*boom.Radius^2];

% MOI check
    % syms x y z;
    Iint0 = crossr([x,y,z]')'*crossr([x,y,z]');
    % Iint1 = int(Iint0, z, -sqrt(boom.Radius^2 - y^2), sqrt(boom.Radius^2 - y^2));
    % Iint2 = int(Iint1, y, -boom.Radius, boom.Radius);
    % I_b_i = double(int(boom.rho*Iint2, x, 0, boom.Length));

% First moment of mass
% ADDED A NEGATIVE SIGN HERE
c_B_i = boom.m*[boom.Length/2; 0; 0];
c_B_ix = crossr(c_B_i);

% Put these constants into boom structure
boom.c_B_i = c_B_i;
clear c_B_i
boom.c_B_ix = c_B_ix;
clear c_B_ix
boom.I_b_i = I_b_i;
clear I_b_i
boom.G_B_i = G_B_i;
clear G_B_i
boom.H_B_i = H_B_i;
clear H_B_i
boom.P_B_i = P_B_i;
clear P_B_i

% States ---------------------------------------------------------------- %
% Define individual states
r_B_ia = state(3, 'Position', [0 0 0]');
v_B_ia = state(3, 'Velocity', [0 0 0]');
p_ia = state(9, 'Vectorized DCM',[1 0 0 0 1 0 0 0 1]' );
w_ia = state(3, 'Angular Velocity', [0 0 0]');
%q_eps = state(2*n_psi, 'elastic_coordinates');
%q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [1.0 0 0 0.0 0 0]');
q_eps_dot = state(2*n_psi, 'elastic_coordinates', [0 0 0 0 0 0]' );


initial_condition_case = 1;
% 1 - Step response zero initial conditions
% 2 - one boom deflected
% 3 - opposition booms deflected
% 4 - adjacent booms deflected
% 5 - Free response deflection in each direction for every boom
% REMEMBER TO CHANGE TENSIONS AS WELL

q0_deform = [-1.43490183425856e-23 , -3.07568735419596e-24, -6.55515369732691e-24, 1.98675611333128, -1.11896251484622, 0.313876452375579]';
switch initial_condition_case
    case 1
        switch comp_number
            case 2
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 3
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 4
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 5
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            otherwise
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
        end
        
        
        % steady state boom actuation initial conditions
    case 2
        switch comp_number
            case 2
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, q0_deform);
            case 3
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 4
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 5
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            otherwise
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
        end
    case 3
        % steady state boom actuation initial conditions
        switch comp_number
            case 2
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, q0_deform);
            case 3
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 4
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, q0_deform);
            case 5
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            otherwise
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
        end
    case 4
        % steady state boom actuation initial conditions
        switch comp_number
            case 2
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, q0_deform);
            case 3
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, q0_deform);
            case 4
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            case 5
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
            otherwise
                q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [0 0 0 0 0 0]');
        end
    case 5
        % Free response
        % UPDATE THE TENSIONS
        q_eps = state_2D_elastic_coordinate(6, "elastic coordinates", boom.PSI, boom.Length, [1 0 0 1 0 0]');
        
    otherwise
end

% Assemble componet state vector and velocities
q_B1 = state_vector({r_B_ia, p_ia, q_eps});
vel_B1 = state_vector({v_B_ia, w_ia, q_eps_dot});

% Create component object ----------------------------------------------- %
comp = component(comp_number, q_B1, vel_B1, @build_M, @build_K_bar, @build_fnon, @build_f, @y2qv_and_q_dot, boom, @build_Gamma);
if nargin == 0
    tf = 60 ;
    dtN = tf*10;
    comp = comp.simulate_component(tf, dtN);
    %comp = comp.simulate_component(tf, 10);
    comp.state.plot;
    comp.compute_energies()
    %comp.velocity.plot;
end

end







