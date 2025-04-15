classdef state_vector

    properties
        states
        number_of_states
        dimensions
        time
        history
        y0
    end % end properties

    methods
        % Constructor
        function obj = state_vector(states)
            if (nargin > 0)
                obj.states = states;
                obj.number_of_states = length(states);
                obj.dimensions = zeros(1,obj.number_of_states);
                obj.y0 = [];
                
                for lv1 = 1: obj.number_of_states
                    obj.dimensions(1,lv1) = states{lv1}.dimension;
                end

                if ( exist('states{1}.y0', 'var')|| ~isempty(states{1}.y0) )
                    for lv1 = 1: obj.number_of_states
                        obj.y0 = [obj.y0; states{lv1}.y0];
                    end
                end
                
            end
        end

        % Set initial condition of state vector
        function obj = set_initial_condition(obj, y0)
            % Set overall initial condition
            obj.y0 = y0;

            % Decompose initial conditions
            decomp = obj.decompose(y0);

            % Set individual initial conditions
            for lv1 = 1:obj.number_of_states
                obj.states{lv1}.y0 = decomp{lv1};
            end
        end
        
        % Import state history and parse through it
        function obj = import_history(obj, time, history)
                obj.time =time;
                obj.history = history;

                i = 0;
                for lv1 = 1:obj.number_of_states
                    histroy_i = history(i+1:obj.dimensions(lv1)+i, :);
                    obj.states{lv1} = obj.states{lv1}.import_history(time,histroy_i);
                    i = i + obj.dimensions(lv1);
                end
        end
        
        % Output a cell array of the states
        function decomp = decompose(obj,q)
            decomp = cell(1, obj.number_of_states);
            n_q = 1;
            for lv1 = 1:obj.number_of_states
                decomp{lv1} = q(n_q: sum(obj.dimensions(1:lv1)) );
                n_q = n_q + obj.dimensions(lv1);
            end
        end
        
        % Call respective plot functions
        function plot(obj)
            for lv1 = 1:obj.number_of_states
                obj.states{lv1}.plot();
            end

        end
    end % end methods
end % end class














