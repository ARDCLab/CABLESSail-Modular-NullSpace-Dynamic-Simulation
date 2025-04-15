classdef (Abstract) component_assembly_BASE
    properties
        % Default value if not overwritten by user
        %ode_options = odeset('RelTol',1e-6,'AbsTol',1e-7);
        %ode_options = odeset('RelTol',1e-9,'AbsTol',1e-12);
        ode_options = odeset('RelTol',1e-12,'AbsTol',1e-13);

        % Cell array of assembly components ----------------------------
        components

        % Coodriantes --------------------------------------------------
        coordinate
        velocity
        velocity_reduced

        % Energies -----------------------------------------------------
        PE
        KE
        TE

        % simulation outputs -------------------------------------------
        t 
        y 
        
    end % end properties

    methods (Abstract)
        Upsilon     = compute_Upsilon(obj, t, coords, vels_hat, q, v_hat)  
        Upsilon_dot = compute_Upsilon_dot(obj, t, coords, vels_hat, q, v_hat)  
        [coords, q] = q_hat2q(obj, coords_reduced, q_hat)
        f_gen_ext = compute_f_gen_ext(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array);
    end

    methods
        % Constructor
        function obj = component_assembly_BASE(components, velocity_reduced, varargin )
            %component_assembly_BASE Subclass constructor
            %   Assign the component array and reduced velocty vector to
            %   the class property.
            %
            % --------------------------------------------------------- %
            
            % I think this constructuor can be used by every assembly...
            
            % Assign the input component list to the component property
            obj.components = components;
            obj.velocity_reduced = velocity_reduced;

            if nargin >2
                % derp. Profit?
            end

            % Assemble the velocity and coordinates
            coordinate = {};
            velocity = {};

            for i = 1:length(components)
                coordinate = [ coordinate(:)', obj.components{i}.coordinate.states(:)' ];
                velocity = [ velocity(:)', obj.components{i}.velocity.states(:)' ];
            end
            obj.coordinate = state_vector(coordinate);
            obj.velocity = state_vector(velocity);

        end
            
        dydt = assembly_ode(obj, t, y) % ode function for simulation
        obj  = set_initial_conditions(obj, t, coords_reduced_0, q_hat_0, vels_reduced_0, v_hat_0)


    end % end methods

    methods
        function [coords, vels, vels_hat, q, v, v_hat, Upsilon, Upsilon_dot] = decomepose_state(obj, y)
            %decomepose_state Decompose a current total state column vector
            %into the individual coordinates and velocities.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Get the length of total coordiante and velocity vector
            n_coord = sum(obj.coordinate.dimensions);
            n_vel_reduced = sum(obj.velocity_reduced.dimensions);

            % Check that the lengths match. If not, throw an error
            if length(y) ~= (n_coord+n_vel_reduced)
                error('Component state decomposition failed. Length of input does not match the total state length of the component.')
            end
            
            % Compute the state vectors (doubles)
            q           = y(1:n_coord);
            v_hat       = y(n_coord+1:end);

            % Decompose the vectors into cell arrays of individual states
            coords      = obj.coordinate.decompose( q );
            vels_hat    = obj.velocity_reduced.decompose( v_hat );

            % Compute velocity state vector
            Upsilon     = obj.compute_Upsilon(0, coords, vels_hat, q, v_hat);
            Upsilon_dot = obj.compute_Upsilon_dot(0, coords, vels_hat, q, v_hat);
            v           = Upsilon*v_hat;
            vels        = obj.velocity.decompose( v );

        end

        function [comp_coords_array, comp_vels_array, comp_q_array, comp_v_array] = decompoese_component_states(obj, coords, vels, q, v)
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Determine the number of components for array size
            number_of_components = length(obj.components);
            
            % variables to keep track variable sizes
            ncoords    = 0;
            nvels      = 0;
            nq         = 0;
            nv         = 0;

            % Create empy arrays for the variables
            comp_coords_array   = cell(1,number_of_components);
            comp_vels_array     = cell(1,number_of_components);
            comp_q_array        = cell(1,number_of_components);
            comp_v_array        = cell(1,number_of_components);
            comp_y_array        = cell(1,number_of_components);

            % Loop through components
            for lv1 = 1:number_of_components
                % Get current component
                comp_i = obj.components{lv1};
                
                % Get variable lengths
                nq_i        = sum(comp_i.coordinate.dimensions);
                nv_i        = sum(comp_i.velocity.dimensions);
                ncoords_i   = comp_i.coordinate.number_of_states;
                nvels_i     = comp_i.velocity.number_of_states;

                % Fill in variable arrys
                comp_coords_array_i = coords(1+ncoords : ncoords+ncoords_i);
                comp_vels_array_i   = vels(1+nvels : nvels+nvels_i);
                comp_q_array_i      = q(1+nq : nq+nq_i);
                comp_v_array_i      = v(1+nv : nv+nv_i);
                comp_y_array_i      = [comp_q_array_i; comp_v_array_i];

                % Add to arrays
                comp_coords_array{lv1}   = comp_coords_array_i;
                comp_vels_array{lv1}     = comp_vels_array_i;
                comp_q_array{lv1}        = comp_q_array_i;
                comp_v_array{lv1}        = comp_v_array_i;
                comp_y_array{lv1}        = comp_y_array_i;

                % Update lengths
                nq      = nq + nq_i;
                nv      = nv + nv_i;
                nvels   = nvels + nvels_i;
                ncoords = ncoords + ncoords_i;
            end % end for
        end % end function

        function obj = simulate(obj, tf, dtN)
            %simulate Simulates the component with the initial
            %conditions specific in the state properties.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            tSteps = linspace(0,tf,dtN);
            y_0 = [obj.coordinate.y0;
                   obj.velocity_reduced.y0];
            [tout, yout] = ode15s(@(t,y)obj.assembly_ode(t, y) , tSteps', y_0, obj.ode_options);
            obj.t = tout';
            obj.y = yout';

            obj.coordinate  = obj.coordinate.import_history(tout', yout(:, 1:sum(obj.coordinate.dimensions))' );
            obj.velocity_reduced    = obj.velocity_reduced.import_history(tout', yout(:, sum(obj.coordinate.dimensions)+1:end)' );
        end

        function K = compute_K(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array)
            %compute_K Compute the aggregated stiffness matrix of the total
            %system for the current state.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Build the stiffness matrix of the first body
            K1 = obj.components{1}.compute_K(t, comp_coords_array{1}, comp_vels_array{1}, comp_q_array{1}, comp_v_array{1});

            % Loop through all remaining components and assembly K matrix
            K = K1;
            for lv1 = 2:length(obj.components)
                Ki = obj.components{lv1}.compute_K(t, comp_coords_array{lv1}, comp_vels_array{lv1}, comp_q_array{lv1}, comp_v_array{lv1});
                K = blkdiag(K, Ki);
            end
        end

        function M = compute_M(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array)
            %Compute_M_bar Compute the aggregated mass matrix of the total
            %system for the current state.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Build the stiffness matrix of the first body
            M1 = obj.components{1}.compute_M(t, comp_coords_array{1}, comp_vels_array{1}, comp_q_array{1}, comp_v_array{1});

            % Loop through all remaining components and assembly K matrix
            M = M1;
            for lv1 = 2:length(obj.components)
                Mi = obj.components{lv1}.compute_M(t, comp_coords_array{lv1}, comp_vels_array{lv1}, comp_q_array{lv1}, comp_v_array{lv1});
                M = blkdiag(M, Mi);
            end
        end

        function f_non = compute_f_non(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array)
            %compute_f_non Compute the aggregated non-linear function
            %vector of the total system for the current state.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Build the stiffness matrix of the first body
            f_non1 = obj.components{1}.compute_f_non(t, comp_coords_array{1}, comp_vels_array{1}, comp_q_array{1}, comp_v_array{1});

            % Loop through all remaining components and assembly K matrix
            f_non = f_non1;
            for lv1 = 2:length(obj.components)
                f_non_i = obj.components{lv1}.compute_f_non(t, comp_coords_array{lv1}, comp_vels_array{lv1}, comp_q_array{lv1}, comp_v_array{lv1});
                f_non = [f_non; f_non_i];
            end
        end

        function f_gen = compute_f_gen(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array) 
            %compute_f Compute the aggregated generalized force vector of
            %the total system for the current state.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Build the stiffness matrix of the first body
            f_gen1 = obj.components{1}.compute_f_gen(t, comp_coords_array{1}, comp_vels_array{1}, comp_q_array{1}, comp_v_array{1});

            % Loop through all remaining components and assembly K matrix
            f_gen = f_gen1;
            for lv1 = 2:length(obj.components)
                f_gen_i = obj.components{lv1}.compute_f_gen(t, comp_coords_array{lv1}, comp_vels_array{lv1}, comp_q_array{lv1}, comp_v_array{lv1});
                f_gen = [f_gen; f_gen_i];
            end
        end

        function Gamma = compute_Gamma(obj, t, comp_coords_array, comp_vels_array, comp_q_array, comp_v_array)
            %compute_Gamma Compute the aggregated gamma matrix of the entire
            %system for the current state.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %
            
            % Build the stiffness matrix of the first body
            Gamma1 = obj.components{1}.compute_gamma(t, comp_coords_array{1}, comp_vels_array{1}, comp_q_array{1}, comp_v_array{1});

            % Loop through all remaining components and assembly K matrix
            Gamma = Gamma1;
            for lv1 = 2:length(obj.components)
                Gamma_i = obj.components{lv1}.compute_gamma(t, comp_coords_array{lv1}, comp_vels_array{lv1}, comp_q_array{lv1}, comp_v_array{lv1});
                Gamma = blkdiag(Gamma, Gamma_i);
            end
        end
        
        function obj = import_coordinates_velocities(tout, yout)
            %import_coordinates Assign time history values to the system
            %coordiantes and velocities, and assign the respect time
            %history to the compone t coordinates and velocities.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % Assign histroy to coordinates and velocity

            % Loop through components, assign coordinate and velocity
            % histories to each resepctive component. Use the component
            % import function

        end

        function obj = compute_energies()
            %compute_energies Compute the potential, kinetic, and total
            %energy of the total system from the state history based on the
            %mass matrix and stiffness matrix energy calculations.
            %   Detailed explanation goes here
            %
            % ----------------------------------------------------------- %

            % -------------------------- %
            % LOW PRIORITY UNTIL RUNNING %
            % -------------------------- %

            % Compute energy two ways. One with aggregated mass and
            % stiffness matrices, and the other by using the compute
            % energy function of each component
        end

    end
end % end class