classdef bus_component < component_BASE

    properties
        coordinate = state_vector({state(3, 'Position',[0 0 0]'), ...
                        state(9, 'Vectorized DCM', [1 0 0 0 1 0 0 0 1]')});
        velocity = state_vector({state(3, 'Velocity', [0 0 0]'), ...
                        state(3, 'Angular Velocity', [0 0 0]')});
    end % properites

    methods

        function obj = bus_component()
            % Call the super class constructor
            obj = obj@component_BASE();
            obj = build_constants_structure(obj);
        end

        function obj = build_constants_structure(obj)
            % Dimensions of the bus
            const.width = 0.30; % meter
            const.depth = 0.30; % meter
            const.height = 0.10; % meter

            % Mass properties
            const.mass_bus = 1; % kg
            const.mass_sail = 50; % kg
            const.mass_total = const.mass_bus + const.mass_sail;
            
            % Inertia properties
            a = 29.5*sqrt(2);
            I_sail = 1/12*const.mass_sail*diag([a^2, a^2, 2*a^2]);
            I_bus = 1/12*const.mass_bus*diag([const.height^2+const.depth^2, const.depth^2+const.width^2, const.height^2+const.width^2]);
            const.I_Bbb_b = I_sail + I_bus;

            % Assign the constant to the structure
            obj.const = const;
        end


    end


end % classdef