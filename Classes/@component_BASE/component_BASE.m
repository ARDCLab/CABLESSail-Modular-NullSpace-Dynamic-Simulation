classdef (Abstract) component_BASE
    
    properties
        ode_options = odeset('RelTol',1e-3,'AbsTol',1e-6); 
        % Default value if not overwritten by user
    end

    properties (Abstract) % MUST be overwritten by subclass
        coordinate  % Cell arrays of state object coordinates
        velocity    % Cell arrays of state object velocities
    end

    properties
        % Component constant structure -----------------------------------
        const 

        % Energies -------------------------------------------------------
        KE  % Row vector of kinetic energy values
        PE  % Row vector of potential energy values
        TE  % Row vector of total energy values

        % simulation outputs ---------------------------------------------
        t   % Row vector of time values
        y   % Matrix of column vectors containing state values (not state objects)

    end % properites

    methods (Abstract) % MUST be defined by subclass
        M = compute_M(obj, t, coords, vels, q, v)         % Mass matrix
        K = compute_K(obj, t, coords, vels, q, v)         % Stiffness matrix
        f_non = compute_f_non(obj, t, coords, vels, q, v) % Nonlinear function
        f_gen = compute_f_gen(obj, t, coords, vels, q, v) % Generalized forces
        gamma = compute_gamma(obj, t, coords, vels, q, v) % Attitude rate mapping matrix. Add default that is just identity?
    end

    methods
        function [coords, vels, q, v] = decomepose_state(obj, y)
            %decomepose_state Decompose a current total state column vector
            %into the individual coordinates and velocities.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Get the length of total coordiante and velocity vector
            n_coord = sum(obj.coordinate.dimensions);
            n_vel = sum(obj.velocity.dimensions);

            % Check that the lengths match. If not, throw an error
            if length(y) ~= (n_coord+n_vel)
                error('Component state decomposition failed. Length of input does not match the total state length of the component.')
            end
            
            q = y(1:n_coord);
            v = y(n_coord+1:end);
            coords = obj.coordinate.decompose( q );
            vels = obj.velocity.decompose( v );
        end


        function obj = component_BASE()
            %component_BASE Base constructor.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % IDK what needs to go here for the base class that should
            % never be called....
        end

        function dydt = component_ode(obj, t, y)
            %component_ode Compute dy/dt for the total system dynamics.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %
            % Extract coordinates and velocities from the current state
            [coords, vels, q, v] = obj.decomepose_state( y );
            
            % Compute dynamic matrices
            % TODO: pass in the outputs of decompoese_state so it doesn't
            % have to be called inside every compute function ~for speeed~
            f_non   = obj.compute_f_non(t, coords, vels, q, v);
            gamma   = obj.compute_gamma(t, coords, vels, q, v);
            K       = obj.compute_K(t, coords, vels, q, v);
            M       = obj.compute_M(t, coords, vels, q, v);

            % Compute M*v_dot = ()
            M_times_v_dot = -f_non - gamma'*K*q;
            
            % Multiply through my inverse of M and solve for v_dot
            v_dot = pinv(M) * M_times_v_dot;

            % Check Natural Frequencies --------------------------------- %
                % MinvK = obj.gamma(t,y,obj.const)*pinv(obj.M_bar(t,y,obj.const))*obj.gamma(t,y,obj.const)'*obj.K(t,y,obj.const);
                % w_nat = sqrt(eig(MinvK));
                % w_nat_q = sqrt(eig(MinvK(end-5:end, end-5:end)));
            % ----------------------------------------------------------- %

            % Return dy/dt
            q_dot = gamma*v;
            dydt = [q_dot; v_dot];
        end

        function obj = simulate_component(obj, tf, dtN)
            %simulate_component Simulate the full system dynamics with the
            %initial conditions from the coordinates and velocities.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            tSteps = linspace(0,tf,dtN);
            
            y_0 = [obj.coordinate.y0;
                   obj.velocity.y0];

            %[tout, yout] = ode45(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
            [tout, yout] = ode113(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
            
            % Record output to the object
            obj.t = tout';
            obj.y = yout';
            
            % Import history to state and velocity objects
            obj.coordinate  = obj.coordinate.import_history(tout', yout(:, 1:sum(obj.coordinate.dimensions))' );
            obj.velocity    = obj.velocity.import_history(tout', yout(:, sum(obj.coordinate.dimensions)+1:end)' );

            % TODO: Create set or import methods for coordinates and
            % velocities that then copy them to the respective coordinates
            % and velocities of the components
        end

        function obj = import_coordinates_velocities(tout, yout)
            %import_coordinates Assign time history values to the component
            %coordiante and velocity properties
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %
            % Check that the dimensions matchup

            % Assign histroy to coordinates and velocity assuming they are
            % in the expected order
        end

        function obj = compute_energies(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            obj.PE = zeros(size(obj.t));
            obj.KE = zeros(2, length(obj.t));

            for lv1 = 1:length(obj.t)
                vi = obj.velocity.history(:,lv1);
                xi = obj.coordinate.history(:,lv1);
                yi = [xi; vi];
                
                % comment out for the hub
                % obj.KE(2,lv1) = KE_velocity(xi, vi, obj.const);

                obj.KE(1,lv1) = 1/2 * vi' * obj.M_bar(obj.t(lv1), yi, obj.const )*vi;
                obj.PE(1,lv1) = 1/2 * xi' * obj.K(obj.t(lv1), yi, obj.const )*xi;

                % if obj.KE(1,lv1) < 0
                %     disp('Warning, negative kinetic energy')
                % end
            end
            obj.TE = obj.KE(1,:) + obj.PE;
            obj.TE = abs(obj.TE -obj.TE(1,1));


            figure()
            subplot(3,1,1)
            plot(obj.t, obj.PE)
            grid on
            title('Potential Energy', 'FontSize', 24)
            subplot(3,1,2)
            plot(obj.t, obj.KE(1,:))
            grid on
            hold on
            title('Kinetic Energy', 'FontSize', 24)
            subplot(3,1,3)
            plot(obj.t, obj.TE)
            grid on
            title('Change in Total Energy', 'FontSize', 24)
        end
    end % methods
end % classdef