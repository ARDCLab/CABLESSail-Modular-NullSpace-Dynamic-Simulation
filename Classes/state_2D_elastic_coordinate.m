classdef state_2D_elastic_coordinate < state

    properties
        % INHERITED:
        % dimension
        % description
        % time
        % history

        % NEW:
        basis_functions
        deflection % input x, output [u1(x), u2(x)]' for all t
        Length

    end % end properties

    methods
        % The constructor is inherited????
        % Constructor
        % function obj = state_2D_elastic_coordinate(dim, descrp, basis_functs, L, IC)
        function obj = state_2D_elastic_coordinate(varargin)
            obj@state(varargin{:,:});
        end

        function obj = define_basis(obj, basis_functs, L)
            obj.basis_functions = basis_functs;
            obj.Length = L;
        end

        function obj = import_history(obj, time, history)
            obj = import_history@state(obj, time, history);
            obj.deflection = obj.basis_functions * history;
            obj.deflection = matlabFunction(obj.deflection);
            if sum(history)==0
                obj.deflection = @(x) x*zeros(2,length(time));
            end
        end
        

        % OVERRRIDE PLOT?
        function plot(obj)
            % figure()
            % for lv1 = 1:obj.dimension
            %     subplot(3,2,lv1)
            %     plot(obj.time , obj.history(lv1, :))
            %     hold on 
            %     grid on
            % end
            % sgtitle('Elastic Coordinates', 'interpreter','Latex', 'FontSize', obj.TitleSize)
            
            figure()
            % Breaks here if the deflection is ALL 0
            tip_deflection = obj.deflection(obj.Length);
            subplot(2,4,1:2)
            plot(obj.time, tip_deflection(1,:))
            grid on
            ylabel('$u_1$ (m)', 'interpreter','Latex', 'FontSize',obj.LabelSize)

            subplot(2,4,5:6)
            plot(obj.time, tip_deflection(2,:))
            grid on
            ylabel('$u_2$ (m)','interpreter','latex', 'FontSize',obj.LabelSize)
            xlabel('Time (s)', 'interpreter','latex', 'FontSize',obj.LabelSize)
            sgtitle('Tip Deformation', 'interpreter','Latex', 'FontSize', obj.TitleSize)

            subplot(2,4, [3,4,7,8])
            plot(tip_deflection(1,:), tip_deflection(2,:))
            grid on
            ylabel('$u_2$ (m)','interpreter','latex', 'FontSize',obj.LabelSize)
            xlabel('$u_1$ (m)', 'interpreter','Latex', 'FontSize',obj.LabelSize)
            axis equal

            

        end
    end % end methods
end % end class