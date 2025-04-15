classdef state

    properties
        dimension
        description
        time
        history
        TitleSize = 24;
        LabelSize = 18;
        y0
    end % end properties

    methods
        % Constructor
        function s = state(varargin)
            % dim, descp, IC
            if (nargin == 3)
                dim = varargin{1};
                descp = varargin{2};
                IC = varargin{3};
                s.dimension   = dim;
                s.description  = descp;
                if ~isempty(IC)
                    s.y0 = IC;
                end
            end
            if (nargin == 2)
                dim = varargin{1};
                descp = varargin{2};
                s.dimension   = dim;
                s.description  = descp;
            end
        end

        function obj = set_initial_condition(y0)
            obj.y0 = y0;
        end
        
        % Import state history and parse through it
        function obj = import_history(obj, time, history)
            obj.time = time;
            obj.history = history;
        end


        function plot(obj)
            figure()

            if obj.dimension == 3
                for lv1 = 1:obj.dimension
                    subplot(obj.dimension,1,lv1)
                    plot(obj.time , obj.history(lv1, :))
                    hold on
                    grid on
                end
                sgtitle(obj.description, 'FontSize', obj.TitleSize)
                xlabel('Time (s)','FontSize',obj.LabelSize)
            elseif obj.dimension == 6
                for lv1 = 1:obj.dimension
                    subplot(3,2,lv1)
                    plot(obj.time , obj.history(lv1, :))
                    hold on 
                    grid on
                end
                sgtitle('Elastic Coordinates', 'FontSize', obj.TitleSize)
            elseif obj.dimension == 9
                for lv1 = 1:obj.dimension
                    subplot(4,3,lv1)
                    plot(obj.time , obj.history(lv1, :))
                    hold on
                    grid on
                end
                for lv1 = 1:length(obj.time)
                    % [C_ab, ~] = build_C_ab(obj.history(:, lv1), false);
                    % Cba = C_ab';
                    p_ia = obj.history(:, lv1);
                    Cia = p2DCM(p_ia);

                    subplot(4,3,10:12)
                    plot(obj.time(lv1), abs(1-det(Cia)),'b.' )
                    hold on
                    grid on
                end
                title('Determinant-1', 'FontSize', obj.TitleSize)
                sgtitle('DCM Entries', 'FontSize', obj.TitleSize)

                figure()
                hold on
                for lv1 = 1:length(obj.time)
                    % [C_ab, ~] = build_C_ab(obj.history(:, lv1), false);
                    % Cba = C_ab';
                    rad2deg = 180/pi;
                    p_ia = obj.history(:, lv1);
                    Cia = p2DCM(p_ia);
                    [phii,thetai,psii] = DCM2Euler321(Cia);

                    subplot(3,1,1)
                    plot(obj.time(lv1), phii*rad2deg,'b.')
                    hold on

                    subplot(3,1,2)
                    plot(obj.time(lv1), thetai*rad2deg,'b.')
                    hold on

                    subplot(3,1,3)
                    plot(obj.time(lv1), psii*rad2deg,'b.')
                    hold on

                end
                subplot(3,1,1)
                hold on
                grid on
                ylabel("Roll - Degrees", Interpreter="latex")
                xlim([0, max(obj.time)])

                subplot(3,1,2)
                hold on
                grid on
                ylabel("Yaw - Degrees", Interpreter="latex")
                xlim([0, max(obj.time)])

                subplot(3,1,3)
                hold on
                grid on
                ylabel("Pitch - Degrees", Interpreter="latex")
                xlabel("Time - Seconds", Interpreter="latex")
                xlim([0, max(obj.time)])

                sgtitle('Euler Angles - Radians', Interpreter="latex")
            end

        end
    end % end methods
end % end class














