classdef AMT_component < component_BASE

    properties
        coordinate = state_vector({state(3, 'Position',[0 0 0]'), ...
                        state(9, 'Vectorized DCM', [1 0 0 0 1 0 0 0 1]')});
        velocity = state_vector({state(3, 'Velocity', [0 0 0]'), ...
                        state(3, 'Angular Velocity', [0 0 0]')});
    end % properites

    methods

        function obj = AMT_component()
            % Call the super class constructor
            obj = obj@component_BASE();
            obj = build_constants_structure(obj);
        end

        function obj = build_constants_structure(obj)
            % Dimensions of the AMT
            const.width = 0.30; % meter
            const.depth = 0.30; % meter
            const.height = 0.90; % meter

            % Mass properties of AMT
            const.mass_AMT = 50; % kg
            
            % Inertia properties of AMT
            I_AMT = 1/12*const.mass_AMT*diag([const.height^2+const.depth^2, const.depth^2+const.width^2, const.height^2+const.width^2]);
            const.I_Bbb_b = I_AMT;

            % Assign the constant to the structure
            obj.const = const;
        end


    end


end % classdef