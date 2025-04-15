classdef hub_component < component_BASE

    properties
        coordinate = state_vector({state(3, 'Position',[0 0 0]'), ...
                        state(9, 'Vectorized DCM', [1 0 0 0 1 0 0 0 1]')});
        velocity = state_vector({state(3, 'Velocity', [0 0 0]'), ...
                        state(3, 'Angular Velocity', [0 0 0]')});
    end % properites

    methods

        function obj = hub_component()
            % Call the super class constructor
            obj = obj@component_BASE();
            obj = build_constants_structure(obj);
        end

        function obj = build_constants_structure(obj)
            const.a = 1;
            const.m_hub = 50; % kg
            % height
            const.h = 1; % meter
            % radius
            const.R = 1; % meter
            
            m_sail = 50; % kg
            const.m_sail = m_sail;
            a = 29.5*sqrt(2);
            I_sail = m_sail*1/21*diag([a^2, a^2, 2*a^2]);
            const.m = m_sail + const.m_hub;
            
            
            % Moment of Inertia matrix
            const.I_Bbb_b = const.m_hub*[1/12*const.h^2 + 1/4*const.R^2, 0, 0; ...
                          0, 1/12*const.h^2 + 1/4*const.R^2, 0; ...
                          0, 0, 1/2*const.R^2]+I_sail;

            obj.const = const;
        end


    end


end % classdef