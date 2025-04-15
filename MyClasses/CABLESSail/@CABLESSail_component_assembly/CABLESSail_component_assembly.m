classdef CABLESSail_component_assembly < component_assembly_BASE
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % n/a I guess
        const
        sailMembraneTroque
    end

    methods
        [F_gen_r, F_gen_w] = compute_SS_generalized_forces_PLATE(obj, p_ba, q1, q2, q3, q4)
        [Force_sail, Torque_sail] = Compute_Sail_Force_FlatePlate(obj, p_ba, q1, q2, q3, q4)
        
        function obj = CABLESSail_component_assembly(components, velocity_reduced, varargin )
            %CABLESSail_component_assembly Construct an instance of this class
            %   Call the super class constructor.
            obj = obj@component_assembly_BASE(components, velocity_reduced, varargin );
            obj = build_constants_structure(obj); 
        end

        plot_SRPvTime(obj)

        function obj = build_constants_structure(obj)
            % Orientation of the attached booms
            obj.const.C_cb = eye(3);
            obj.const.C_db = [0 1 0;
                             -1 0 0;
                              0 0 1];
            obj.const.C_eb = [-1 0 0;
                              0 -1 0;
                              0 0 1];
            obj.const.C_fb = [0 -1 0;
                              1 0 0;
                              0 0 1];
            
            % Location of boom attachments
            R = obj.components{1}.const.R;
            obj.const.r_cb_b = R*[1;0;0];
            obj.const.r_db_b = R*[0;1;0];
            obj.const.r_eb_b = R*[-1;0;0];
            obj.const.r_fb_b = R*[0;-1;0];

            % Optical properties of the sail (Solar Cruiser parameters from Heaton and
            % Artusio-Glimpse "An Update to the NASA Reference Solal Sail Thrust
            % Model'' 2015
            obj.const.Opt_prop.P        = 4.53914e-6;
            obj.const.Opt_prop.r_tilde  = 0.91;
            obj.const.Opt_prop.s        = 0.94;
            obj.const.Opt_prop.Bf       = 0.79;
            obj.const.Opt_prop.Bb	    = 0.67;
            obj.const.Opt_prop.ef       = 0.025;
            obj.const.Opt_prop.eb       = 0.27;

            
        end

    end
end