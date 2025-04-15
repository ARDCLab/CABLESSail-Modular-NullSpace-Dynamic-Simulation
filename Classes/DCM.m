classdef DCM
    properties (SetAccess = private)
        C_ba
        frame_a
        frame_b
    end
    methods
        % Constructor method
        function obj = DCM(C_ba, frame_a, frame_b)
            obj.C_ba = C_ba;
            obj.frame_a = frame_a;
            obj.frame_b = frame_b;

            if abs(norm(C_ba) - 1) > 2*eps
                warning("Norm of DCM != 1")
            end
        end

        % Get DCM
        function C_ba = get.C_ba(obj)
            C_ba = obj.C_ba;
        end

        % Update DCM
        function obj = update(obj, C_ba_new)
            obj.C_ba = C_ba_new;
            if abs(norm(C_ba_new) - 1) > 2*eps
                warning("Norm of DCM != 1")
            end
        end
        
        % Compute Gamma matrix
        function Gamma = Gamma(obj)
            p1 = obj.C_ba(:,1);
            p2 = obj.C_ba(:,2);
            p3 = obj.C_ba(:,3);
            Gamma = [zeros(3,1), zeros(3,1), p2;
                    p3, zeros(3,1), zeros(3,1);
                    zeros(3,1), p1, zeros(3,1)];
        end

        % Compute S matrix
        function S = S(obj)
            S = (obj.Gamma())';
        end
        
        % Get p vector from DCM
        function p = p(obj)
            p = [obj.C_ba(:,1) ; obj.C_ba(:,2) ; obj.C_ba(:,3)];
        end
    end
end

