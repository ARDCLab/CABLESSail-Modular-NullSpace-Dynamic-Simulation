classdef component_assembly
    properties (Constant)
        ode_options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        %ode_options = odeset('RelTol',1e-7,'AbsTol',1e-9);
        %ode_options = odeset('RelTol',1e-13,'AbsTol',1e-14);
        %ode_options = odeset('RelTol',1e-13,'AbsTol',1e-15);
    end
    properties
        components

        % Boodriantes
        state
        velocity
        y2q_v_qdot
        
        % Null space change of variables matrix 
        % Upsilon
        % Upsilon_dot
        
        % Function handle lists 
        M_bar_list
        K_list
        f_non_list
        f_list
        i_component_list
        Gamma_list
        const_list

        % Energies
        PE
        KE
        TE

        % simulation outputs ---------------------------------------------
        t 
        y 
        

    end % end properties

    methods
        % Constructor
        function cmpnt = component_assembly(components, gen_velocity )
            if (nargin > 0)
                
                cmpnt.components = components; 
                cmpnt.state = state_vector([]);
                cmpnt.velocity = gen_velocity;
                
                % Initialize function handle cell arrays
                cmpnt.M_bar_list        = cell(1, length(components));
                cmpnt.K_list            = cell(1, length(components));
                cmpnt.f_non_list        = cell(1, length(components));
                cmpnt.f_list            = cell(1, length(components)); 
                cmpnt.i_component_list  = cell(1, length(components));
                cmpnt.Gamma_list        = cell(1, length(components));
                cmpnt.const_list        = cell(1, length(components));
                state_list              = [];
                

                % Fill function handle cell arrays
                for lv1 = 1:length(components)
                    cmpnt.M_bar_list{lv1}       = components{lv1}.M_bar;
                    cmpnt.K_list{lv1}           = components{lv1}.K;
                    cmpnt.f_non_list{lv1}       = components{lv1}.f_non;
                    cmpnt.f_list{lv1}           = components{lv1}.f;
                    cmpnt.i_component_list{lv1} = components{lv1}.i_component;
                    cmpnt.Gamma_list{lv1}       = components{lv1}.gamma;
                    cmpnt.const_list{lv1}       = components{lv1}.const;
                    state_list                  = [state_list components{lv1}.state.states];
                end
                cmpnt.state = state_vector(state_list);
            end
            % cmpnt.Upsilon_dot    = Upsilon_dot;
            % cmpnt.Upsilon       = null_space_matrix;

        end


        % NOT IN USE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Extract specific states
        function [r_ia_a, v_ia_a, p_ia, p_ia_dot, omega_ia_i, C_ia, q_eps, q_eps_dot ] = extract_state_components(obj,q,nu)
                r_ia_a = q(1:3,1);
                v_ia_a = nu(1:3,1);
                p_ia = q(4:12,1);
                q_eps = q(13:end,1);
                q_eps_dot = nu(7:end);
                omega_ia_i = nu(4:6,1);
                p_ia_dot = build_p_dot(p_ia, omega_ia_i);
                [C_ai, p_ia] = build_C_ab(p_ia, false);
                C_ia = C_ai';
        end

        % Simulate
        function obj = simulate_component(obj, tf, dtN)
            tSteps = linspace(0,tf,dtN);

            y_0 = [obj.state.y0;
                   obj.velocity.y0];

            [tout, yout] = ode15s(@(t,y)obj.assembly_ode(t, y) , tSteps', y_0, obj.ode_options);
            
            obj.t = tout';
            obj.y = yout';

            obj.state    = obj.state.import_history(tout', yout(:, 1:sum(obj.state.dimensions))' );
            obj.velocity = obj.velocity.import_history(tout', yout(:, sum(obj.state.dimensions)+1:end)' );
        end

    end % end methods

    % STATIC METHODS
    % TODO: Make these constant methods that are a part of the default
    % class that are never overwritten by the user
    % ------------------------------------------------------------------- %
    methods
        function K = compute_K(obj, t, y_cell)
            yi = y_cell{1};
            K = obj.K_list{1}(t, yi, obj.const_list{1}) ;
            for lv1 = 2:length(obj.components)
                yi = y_cell{lv1};
                K = blkdiag( K, obj.K_list{lv1}(t, yi, obj.const_list{lv1}) );
            end
        end

        function M_bar = compute_M_bar(obj, t, y_cell)
            yi = y_cell{1};
            M_bar = obj.M_bar_list{1}(t, yi, obj.const_list{1});
            for lv1 = 2:length(obj.components)
                yi = y_cell{lv1};
                M_bar = blkdiag( M_bar, obj.M_bar_list{lv1}(t, yi, obj.const_list{lv1}) );
            end
        end

        function f_non = compute_f_non(obj, t, y_cell)
            yi = y_cell{1};
            f_non = obj.f_non_list{1}(obj, t, yi);
            for lv1 = 2:length(obj.components)
                yi = y_cell{lv1};
                f_non = [ f_non; obj.f_non_list{lv1}(t, yi, obj.const_list{lv1}) ];
            end
        end

        function f = compute_f(obj, t, y_cell)           
            yi = y_cell{1};
            f = obj.f_list{1}(obj, t, yi);
            for lv1 = 2:length(obj.components)
                yi = y_cell{lv1};
                f = [f; obj.f_list{lv1}(t, yi, obj.const_list{lv1})];
            end
        end

        function Gamma = compute_Gamma(obj, t, y_cell)
            yi = y_cell{1};
            Gamma = obj.Gamma_list{1}(t, yi, obj.const_list{1});
            for lv1 = 2:length(obj.components)
                yi = y_cell{lv1};
                Gamma = blkdiag( Gamma, obj.Gamma_list{lv1}(t, yi, obj.const_list{lv1}) );
            end
        end


    end
end % end class