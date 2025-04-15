classdef null_space_system
    properties
        components %an array of compoent objects
        constrained_state_dimension % state dimension post null space
        Upsilon % Null space constraining matrix
    end

    methods
        function assemble_unconstrained_system()
        end

        function assemble_constrained_system
        end

        function set_initial_conditions
        end

        function simulate
        end

        function post_process
        end

        function plot
        end

    end

    % is this all worth it? Should I cut back on what this object is doing
    % and just have it assemble the dynamics from the components, then pass
    % off a function handle to the normal ODE solver. It might be easier to
    % do than and THEN turn it into more object oriented programming.
end