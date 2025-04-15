classdef phys_vect

    properties
        vec (3,1) double
        description string
        plot_name string
    end % end properties

    methods
        % Define constructor
        function obj = phys_vect(vec, desc, plot_name)
            obj.vec = vec;
            obj.description = desc;
            obj.plot_name = plot_name;
        end

        % Get skew symmetric matrix
        function rx = r_skew(obj)
            rx = [0,            -obj.vec(1), obj.vec(2);
                  obj.vec(1),   0,          -obj.vec(3);
                 -obj.vec(2),   obj.vec(1), 0];  
        end

        % Update vector
        function obj = update_vec(obj, new_vec)
            obj.vec = new_vec;
        end

    end % end methods
    
end % end class
