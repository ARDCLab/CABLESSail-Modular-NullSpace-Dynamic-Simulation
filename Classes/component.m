classdef component
    
    properties (Constant)
        ode_options = odeset('RelTol',1e-7,'AbsTol',1e-9);
        % ode_options = odeset('RelTol',1e-9,'AbsTol',1e-11);
        %ode_options = odeset('RelTol',1e-13,'AbsTol',1e-15);
    end
    properties
        % size of state, make sure it agrees
     %   state_dimension

        % State and velocity vectors -------------------------------------
        state
        velocity 

        % Energies -------------------------------------------------------
        KE
        PE
        TE

        % Component constant structure -----------------------------------
        const 

        % Component number -----------------------------------------------
        i_component 

        % EoM function handles, depend on the current state, input, and  constants
        M_bar 
        K 
        f_non 
        f 
        y2q_v_qdot 
        gamma 

        % post process function handle -----------------------------------
        post_process 

        % simulation outputs ---------------------------------------------
        t 
        y 

    end % properites
    methods
    % constructor
        function obj = component(i_component, state, velocity, M_bar, K, f_non, f, y2q_v_qdot, const, gamma)
        if nargin > 0    

            % optional gamma matrix
            if nargin == 10
                obj.gamma = gamma;
            else
                obj.gamma = @(t,y,const) 1;
            end

            % Component number
            obj.i_component = i_component;

            % state vectors
            obj.state       = state;
            obj.velocity    = velocity;

            % function handles, depend on the current state and constants
            obj.M_bar       = M_bar;
            obj.K           = K;
            obj.f_non       = f_non;
            obj.f           = f;
            obj.y2q_v_qdot  = y2q_v_qdot;
            obj.const       = const;
        end 
        end

   % ODE of components
        function dydt = component_ode(obj, t, y)
            [q, nu_hat, q_dot] = obj.y2q_v_qdot(t,y,obj.const);
            
            % normalize the DCM --------------------- %
                % Ci = [q(4:6,1), q(7:9,1), q(10:12,1) ];
                % Ci = normalizeDCM(Ci  );
                % q(4:6,1) = Ci(:,1);
                % q(7:9,1) = Ci(:,2);
                % q(10:12,1) = Ci(:,3);
            % --------------------------------------- %

            % Add Artificial Damping ---------------- %
            %delta = 1e-3;
             delta = 0;
            %D = delta*blkdiag(eye(6), eye(6));
            D = zeros(length(nu_hat));
            % --------------------------------------- %

            % if t>10
            %     disp(' ');
            % end

            M_times_v_dot = -obj.f_non(t,y,obj.const) - obj.gamma(t,y,obj.const)'*obj.K(t,y,obj.const)*q;
            M_times_v_dot = M_times_v_dot - D*nu_hat;
            v_dot = pinv(obj.M_bar(t,y,obj.const)) * M_times_v_dot;


            % Check Natural Frequencies --------------------------------- %
                MinvK = obj.gamma(t,y,obj.const)*pinv(obj.M_bar(t,y,obj.const))*obj.gamma(t,y,obj.const)'*obj.K(t,y,obj.const);
                w_nat = sqrt(eig(MinvK));
                w_nat_q = sqrt(eig(MinvK(end-5:end, end-5:end)));
            % ----------------------------------------------------------- %

            dydt = [q_dot;v_dot];
        end
        
    % Simulate the component
        function obj = simulate_component(obj, tf, dtN)
            tSteps = linspace(0,tf,dtN);
            
            y_0 = [obj.state.y0;
                   obj.velocity.y0];

            %[tout, yout] = ode45(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
            %[tout, yout] = ode15s(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
            % [tout, yout] = ode15s(@(t,y)obj.component_ode(t, y) , tSteps', y_0);
            [tout, yout] = ode113(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
            
            % Record output to the object
            obj.t = tout';
            obj.y = yout';
            
            % Import history to state and velocity objects
            obj.state    = obj.state.import_history(tout', yout(:, 1:sum(obj.state.dimensions))' );
            obj.velocity = obj.velocity.import_history(tout', yout(:, sum(obj.state.dimensions)+1:end)' );
        end


    % [DEPRECATED]
    % Simulate the component, looping thorugh dt
        function obj = simulate_component_looped(obj, y_0, tf, dtN)
            tSteps = linspace(0,tf,dtN);
            
            %%%%%%% TEMP TO REMOVE ELASTICS %%%%%%%
            % y_0(13:18) = zeros(6,1);
            % y_0(end-5:end) = zeros(6,1);
            %%%%%%% TEMP TO REMOVE ELASTICS %%%%%%%

            yout = zeros(length(y_0),dtN );
            yout(:,1) = y_0;

            for lv1 = 2:dtN
                %[tout, yout] = ode45(@(t,y)obj.component_ode(t, y) , tSteps', y_0, obj.ode_options);
                [~, yi] = ode15s(@(t,y)obj.component_ode(t, y) , tSteps(lv1-1:lv1), y_0, obj.ode_options);
                yi = yi';

                % normalize the DCM --------------------- %
                Ci = [yi(4:6,1), yi(7:9,1), yi(10:12,1) ];
                Ci = normalizeDCM(Ci);
                yi(4:6,1) = Ci(:,1);
                yi(7:9,1) = Ci(:,2);
                yi(10:12,1) = Ci(:,3);
                % --------------------------------------- %

                yout(:, lv1) = yi(:, end);
            end
            
            tout = tSteps;
            obj.t = tout;
            obj.y = yout;

            obj.state    = obj.state.import_history(tout, yout(1:sum(obj.state.dimensions), :) );
            obj.velocity = obj.velocity.import_history(tout, yout(sum(obj.state.dimensions)+1:end,:) );

        end

    % Compute the energies of the component
        function obj = compute_energies(obj)
            obj.PE = zeros(size(obj.t));
            obj.KE = zeros(2, length(obj.t));

            for lv1 = 1:length(obj.t)
                vi = obj.velocity.history(:,lv1);
                xi = obj.state.history(:,lv1);
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