classdef SailMesh
    %SailMesh Class describing the shape of a solar sail membrane that can
    %be fit to any set of boom tip deformations.
    %
    %   The membrane shape must be defined at declaration with a zero tip
    %   deformation. The XYZ shape is stored within each sail# property,
    %   and is never altered. Copies of this shape must be made whenever
    %   the shape is fit to a set of derofmred booms to compute torques.

    properties (SetAccess = immutable)
        % Immutable means that these can only be set within the constructor  

        % Length of the booms
        L
        % Number of nodes across one dimension of the sail mesh
        n
        % Positions of constraints defining the sail plane
        center = [0 0 0]';
        tip1
        tip2
        tip3
        tip4
        % Structures that contain information describing the shape of the
        % sail mesh for each quadrant of the sail membrane
        sail1
        sail2
        sail3
        sail4
        % CHECK. Not sure how these are used at the moment
        % Old. Used for the flate plate approximation
        plane1
        plane2
        plane3
        plane4
        % Structure containing the optical properties defining the
        % reflection and abosportion properties of the sail
        Opt_prop 
    end
% ----------------------------------------------------------------------- %
%                           CORE FUNCTIONS                                %
% ----------------------------------------------------------------------- %
    methods %-----------------------------------------------------------
        function obj = SailMesh(L, n, MembraneBasisFunctionsCell, functionAmplitudeCell)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            % Creat mesh for all sail quadrants

            obj.Opt_prop = load_optical_properties();

            obj.L = L;
            obj.n = n;

            obj.sail1 = obj.createSailQuadrantMesh(1);
            obj.sail2 = obj.createSailQuadrantMesh(2);
            obj.sail3 = obj.createSailQuadrantMesh(3);
            obj.sail4 = obj.createSailQuadrantMesh(4);

            obj.tip1 = [obj.L 0 0]';
            obj.tip2 = [0 obj.L 0]';
            obj.tip3 = [-obj.L 0 0]';
            obj.tip4 = [0 -obj.L 0]';


            obj.sail1 = applyBasisFunction(obj, 1, MembraneBasisFunctionsCell{1}, functionAmplitudeCell{1});
            obj.sail2 = applyBasisFunction(obj, 2, MembraneBasisFunctionsCell{2}, functionAmplitudeCell{2});
            obj.sail3 = applyBasisFunction(obj, 3, MembraneBasisFunctionsCell{3}, functionAmplitudeCell{3});
            obj.sail4 = applyBasisFunction(obj, 4, MembraneBasisFunctionsCell{4}, functionAmplitudeCell{4});
        end % constructor

        function sailOUT = applyBasisFunction(obj, sailNumber, basisFunctionHandle, functionAmplitude)
            % Apply deflections to each sail quadrant
            % The signs of X and Y are adjusted to be positive so the basis
            % functions behave as they would in the first quadrant.
            switch sailNumber
                case 1
                    sailOUT = obj.sail1;
                    for lv1 = 1:length(obj.sail1.X)
                        sailOUT.Z(lv1) = obj.sail1.Z(lv1) + functionAmplitude * basisFunctionHandle(obj.sail1.X(lv1), obj.sail1.Y(lv1)) ;
                    end
                    % Update elements array
                    sailOUT.elements = updateZelements(obj.sail1.Z, obj.sail1.elements, obj.n);

                case 2
                    sailOUT = obj.sail2;
                    for lv1 = 1:length(obj.sail2.X)
                        sailOUT.Z(lv1) = obj.sail2.Z(lv1) + functionAmplitude * basisFunctionHandle(-obj.sail2.X(lv1), obj.sail2.Y(lv1)) ;
                    end
                    % Update elements array
                    sailOUT.elements = updateZelements(obj.sail2.Z, obj.sail2.elements, obj.n);

                case 3
                    sailOUT = obj.sail3;
                    for lv1 = 1:length(obj.sail3.X)
                        sailOUT.Z(lv1) = obj.sail3.Z(lv1) + functionAmplitude * basisFunctionHandle(-obj.sail3.X(lv1), -obj.sail3.Y(lv1)) ;
                    end
                    % Update elements array
                    sailOUT.elements = updateZelements(obj.sail3.Z, obj.sail3.elements, obj.n);

                case 4
                    sailOUT = obj.sail4;
                    for lv1 = 1:length(obj.sail4.X)
                        sailOUT.Z(lv1) = obj.sail4.Z(lv1) + functionAmplitude * basisFunctionHandle(obj.sail4.X(lv1), -obj.sail4.Y(lv1)) ;
                    end
                    % Update elements array
                    sailOUT.elements = updateZelements(obj.sail4.Z, obj.sail4.elements, obj.n);
            end
        end

        function sail = createSailQuadrantMesh(obj, sailNumber)
        % Creates the mesh of grid points for the resepctive sail quadrant
        % in the form that is compatible with the mesh generator. This mesh
        % is a flat plane and does not take boom tip or membrane deflection
        % into account. This is the nominal mesh.
        %
        % Inputs: 
        %   L = length of booms
        %   n =  number of length segmenets along one edge.
        %
        % ----------------------------------------------------------------------- %
            % Define the linear point density
            dl =  obj.L/obj.n;

            switch sailNumber
                case 1
                    dx = dl;
                    dy = dl;
                case 2
                    dx = -dl;
                    dy = dl;
                case 3
                    dx = -dl;
                    dy = -dl;
                case 4
                    dx = dl;
                    dy = -dl;
            end
            
            % Create the array of points
            elements = cell(1,obj.n+1);
            nY = obj.n+1;
            for lvX = 1:(obj.n+1)
                % Set the size for column of elements
                elements{lvX} = cell(nY,1);
            
                % for lvY = 1:nY 
                %     elements{lvX}{lvY} = [dx*(lvX-1) dy*(lvY-1) 0 ]';
                % end
        
                for lvY = 1:nY 
                    elements{lvX}{lvY} = [(lvX-1) (lvY-1) 0 ]';
                end
            
                nY = nY-1;
            end
            
            % Plot the array of points
            X = [];
            Y = [];
            Z = [];
            for lvX = 1:(obj.n+1)
                for lvY = 1:length(elements{lvX} )
                    X = [X, elements{lvX}{lvY}(1)];
                    Y = [Y, elements{lvX}{lvY}(2)];
                    Z = [Z, elements{lvX}{lvY}(3)];
                end
            end
            
        
            % Create triangle connection matrix
            T = delaunay(X,Y);
            % C = [1, n+1;  n+1, length(X);   length(X), 1;];
            % T2 = delaunayTriangulation(X',Y', C);
            X = X*dx;
            Y = Y*dy;
        
            nY = obj.n+1;
            for lvX = 1:(obj.n+1)
                for lvY = 1:nY 
                    elements{lvX}{lvY} = [dx 0 0; 0 dy 0;0 0 1]* elements{lvX}{lvY};
                end
            
                nY = nY-1;
            end
        
            if length(T) ~= obj.n^2
                error("Incorrect mesh. Redundant mesh triangles created by Delaunay.")
            end

            sail.X = X;
            sail.Y = Y;
            sail.Z = Z;
            sail.T = T;
            sail.elements = elements;    
        end % function

        function [Force_total, Torque_total] = computeSRPForceandTorque(obj, Sun_dir_b, eps1, eps2, eps3, eps4)
            % This function computed the total effect of the SRP on the
            % solar sail by computing the force across the entired
            % discretizied sail, applying the force at the 2D centroid of
            % the discrete sail segement, and replacing all of these force
            % vectors with a total force and couple (or torque) at the
            % center of the sail.
            %
            % Inputs
            %   Sun_dir_b   - Unit vector pointing toward the sun in
            %               the solar sail frame, where the one direction
            %               lies in the direction of the first boom, and 2
            %               direction is along the second boom, such that
            %               the 3rd direction is out of the flat sail plane
            %               towards the nominal sun direction.
            %
            %   eps#        - The out of plane tip deflection of the
            %               respective boom.
            %
            % Outputs
            %   Force_total - 3x1 vector of the total SRP force applied at
            %               the center of the sail
            %
            %   Torque_total - 3x1 vector of the total SRP torque due to
            %               the sail applied at the center of the sail.
            %-------------------------------------------------------------%
            % Apply the tip deformation to the base mesh
            [sail1temp, sail2temp, sail3temp, sail4temp] = obj.applyTipDeflection(eps1, eps2, eps3, eps4);

            % Compute torques contributed from each sail quadrant
            [Force1,Torque1] = MeshSailModel(sail1temp.X, sail1temp.Y, sail1temp.Z, sail1temp.T, Sun_dir_b, obj.Opt_prop);
            [Force2,Torque2] = MeshSailModel(sail2temp.X, sail2temp.Y, sail2temp.Z, sail2temp.T, Sun_dir_b, obj.Opt_prop);
            [Force3,Torque3] = MeshSailModel(sail3temp.X, sail3temp.Y, sail3temp.Z, sail3temp.T, Sun_dir_b, obj.Opt_prop);
            [Force4,Torque4] = MeshSailModel(sail4temp.X, sail4temp.Y, sail4temp.Z, sail4temp.T, Sun_dir_b, obj.Opt_prop);
            Force1 = reshape(Force1,[3,1]);
            Force2 = reshape(Force2,[3,1]);
            Force3 = reshape(Force3,[3,1]);
            Force4 = reshape(Force4,[3,1]);
            Torque1 = reshape(Torque1,[3,1]);
            Torque2 = reshape(Torque2,[3,1]);
            Torque3 = reshape(Torque3,[3,1]);
            Torque4 = reshape(Torque4,[3,1]);

            Torque_total = Torque1 + Torque2 + Torque3 + Torque4;
            Force_total = Force1 + Force2 + Force3 + Force4;
        end % ENF function


        function [sail1, sail2, sail3, sail4] = applyTipDeflection(obj,eps1, eps2, eps3, eps4)
            % Copy the sail variables
            sail1 = obj.sail1;
            sail2 = obj.sail2;
            sail3 = obj.sail3;
            sail4 = obj.sail4;

            tip1 = obj.tip1;
            tip2 = obj.tip2;
            tip3 = obj.tip3;
            tip4 = obj.tip4;
            
            % Update the tips
            tip1(3) = obj.tip1(3) + eps1;
            tip2(3) = obj.tip2(3) + eps2;
            tip3(3) = obj.tip3(3) + eps3;
            tip4(3) = obj.tip4(3) + eps4;

            % Compute the tip positions relative to the old
            tip1Delta = [obj.L 0 eps1]';
            tip2Delta = [0 obj.L eps2]';
            tip3Delta = [-obj.L 0 eps3]';
            tip4Delta = [0 -obj.L eps4]';

            % Create flat plate deflection function for each boom
            plane1 = deflectionBasisFunctionFlatPlate(tip1Delta, tip2Delta);
            plane2 = deflectionBasisFunctionFlatPlate(tip2Delta, tip3Delta);
            plane3 = deflectionBasisFunctionFlatPlate(tip3Delta, tip4Delta);
            plane4 = deflectionBasisFunctionFlatPlate(tip4Delta, tip1Delta);

            % Apply deflections to each sail quadrant
            for lv1 = 1:length(sail1.X)
                % nominal + Plane deflection + billowing + rippling
                sail1.Z(lv1) = sail1.Z(lv1) + plane1(sail1.X(lv1), sail1.Y(lv1));
            
                sail2.Z(lv1) = sail2.Z(lv1) + plane2(sail2.X(lv1), sail2.Y(lv1));
            
                sail3.Z(lv1) = sail3.Z(lv1) + plane3(sail3.X(lv1), sail3.Y(lv1));
            
                sail4.Z(lv1) = sail4.Z(lv1) + plane4(sail4.X(lv1), sail4.Y(lv1));
            end

            % Update elements array
            sail1.elements = updateZelements(sail1.Z, sail1.elements, obj.n);
            sail2.elements = updateZelements(sail2.Z, sail2.elements, obj.n);
            sail3.elements = updateZelements(sail3.Z, sail3.elements, obj.n);
            sail4.elements = updateZelements(sail4.Z, sail4.elements, obj.n);
        end % function
    end %METHODS
    % Methods ----------------------------------------------------------

% ----------------------------------------------------------------------- %
%                           UTILITY FUNCTIONS                             %
% ----------------------------------------------------------------------- %

methods
    function obj = computeSRPthroughAngles(obj)
            % Compare torque and force generated vs flat plate
            Sun_dir_sun_frame = [0;0;1];

            %SIA = 35; % degrees
            %SIA_span = [0:5:25];
            %SIA_span = 25;
            %CLOCK_span = linspace(0, 360, 36);
            
            % Inialize force and torque matrices
            obj.Force_total_plate = zeros(3,length(obj.CLOCK_span),length(obj.SIA_span));
            obj.Torque_total_plate = obj.Force_total_plate;
            
            obj.Force_total_shape =  obj.Force_total_plate;
            obj.Torque_total_shape = obj.Force_total_plate;
            
            for lv_sia = 1: length(obj.SIA_span)
                for lv_clock = 1:length(obj.CLOCK_span)
                    SIA = obj.SIA_span(lv_sia);
                    clock_angle = obj.CLOCK_span(lv_clock);
            
                    C_b_clock = [cosd(clock_angle) -sind(clock_angle) 0;
                                sind(clock_angle) cosd(clock_angle) 0;
                                 0 0 1 ];
                    C_clock_sun = [cosd(SIA) 0 sind(SIA);
                              0 1 0;
                              -sind(SIA) 0 cosd(SIA)];
            
                    C_a_sun = C_b_clock * C_clock_sun;
                    % assuming C_ab = eye(3) at t=0
                    Sun_dir_a = C_a_sun*Sun_dir_sun_frame;
            
                    % Form C_ab from time history
                    C_ab = eye(3);
            
                    % Transform sun into the body frame
                    Sun_dir_b = C_ab' * Sun_dir_a;
            
                    % Bus location
                    Pos_val.hub = [0 0 0]';
            
                    % Call rigid sail model force and torque computation function
                    Pos_val.tip1 = obj.tip1;
                    Pos_val.tip2 = obj.tip2;
                    [Force_t1, Torque_t1] = RigidSailModel(Pos_val,Sun_dir_b,obj.Opt_prop);
            
                    % Call rigid sail model force and torque computation function
                    Pos_val.tip1 = obj.tip2;
                    Pos_val.tip2 = obj.tip3;
                    [Force_t2, Torque_t2] = RigidSailModel(Pos_val,Sun_dir_b,obj.Opt_prop);
            
                    % Call rigid sail model force and torque computation function
                    Pos_val.tip1 = obj.tip3;
                    Pos_val.tip2 = obj.tip4;
                    [Force_t3, Torque_t3] = RigidSailModel(Pos_val,Sun_dir_b,obj.Opt_prop);
            
                    % Call rigid sail model force and torque computation function
                    Pos_val.tip1 = obj.tip4;
                    Pos_val.tip2 = obj.tip1;
                    [Force_t4, Torque_t4] = RigidSailModel(Pos_val,Sun_dir_b,obj.Opt_prop);

                    % Sum force and torque from all sail pieces
                    obj.Torque_total_plate(:,lv_clock, lv_sia) = Torque_t1.total + Torque_t2.total + Torque_t3.total + Torque_t4.total;
                    obj.Force_total_plate(:,lv_clock,lv_sia) = Force_t1.total + Force_t2.total + Force_t3.total + Force_t4.total;
            
                    if ~isreal(obj.Torque_total_plate) || ~isreal(obj.Force_total_plate)
                        disp("Complex Torque or Force")
                    end
                    % Compute for from the complex shape function
                    [Force1,Torque1] = MeshSailModel(obj.sail1.X, obj.sail1.Y, obj.sail1.Z, obj.sail1.T, Sun_dir_b, obj.Opt_prop);
                    [Force2,Torque2] = MeshSailModel(obj.sail2.X, obj.sail2.Y, obj.sail2.Z, obj.sail2.T, Sun_dir_b, obj.Opt_prop);
                    [Force3,Torque3] = MeshSailModel(obj.sail3.X, obj.sail3.Y, obj.sail3.Z, obj.sail3.T, Sun_dir_b, obj.Opt_prop);
                    [Force4,Torque4] = MeshSailModel(obj.sail4.X, obj.sail4.Y, obj.sail4.Z, obj.sail4.T, Sun_dir_b, obj.Opt_prop);
                    obj.Force_total_shape(:,lv_clock, lv_sia) = Force1' + Force2' + Force3' + Force4';
                    obj.Torque_total_shape(:,lv_clock,lv_sia) = Torque1' + Torque2' + Torque3' + Torque4';
                end % END clock angle loop
            end % END SIA loop
        end % FUNCTION

    function surfaceArea = computeSailMeshArea(obj)
                % computeTriMeshArea(sail1.X, sail1.Y, sail1.Z, sail1.T)
                surfaceArea = 0;
                for lv1 = 1:4
                    switch lv1
                        case 1
                            X = obj.sail1.X;
                            Y = obj.sail1.Y;
                            Z = obj.sail1.Z;
                            T = obj.sail1.T;
                        case 2
                            X = obj.sail2.X;
                            Y = obj.sail2.Y;
                            Z = obj.sail2.Z;
                            T = obj.sail2.T;
                        case 3
                            X = obj.sail3.X;
                            Y = obj.sail3.Y;
                            Z = obj.sail3.Z;
                            T = obj.sail3.T;
                        case 4
                            X = obj.sail1.X;
                            Y = obj.sail4.Y;
                            Z = obj.sail4.Z;
                            T = obj.sail4.T;
                    end
                    X = reshape(X, length(X), 1);
                    Y = reshape(Y, length(Y), 1);
                    Z = reshape(Z, length(Z), 1);
                    P = [X, Y, Z];
                    
                    v1 = P(T(:,2), :) - P(T(:,1), :);
                    v2 = P(T(:,3), :) - P(T(:,2), :);
                    
                    cp = 0.5*cross(v1,v2);
                    
                    surfaceArea = surfaceArea + sum(sqrt(dot(cp, cp, 2)));
                end
    end % FUNCTION


end %METHODS
% ----------------------------------------------------------------------- %
%                          PLOTTING FUNCTIONS                             %
% ----------------------------------------------------------------------- %
    % Plotting and animation methods
    methods


        function AnimateSail(obj, t, tip, z_scaling, folder, fileName)
            % Animates the four beam solar sail without sail material. Currently
            % replicating the single current beam 4 times.
            %
            % INPUTS
            % t : 1 x n time vector
            % u : n x 4 array of tip deflections
            % z_scaling : scaling to apply to the z axis of the plottings
            %
            %-------------------------------------------------------------------------
            width = 30;
            height = 20;
            % Initialize video to save animation
            filename = [folder, '/', fileName];
            %myVideo = VideoWriter(filename, 'MPEG-4'); %open video file
            myVideo = VideoWriter(filename, 'MPEG-4'); %open video file
            myVideo.FrameRate = 60;  %can adjust this, 5 - 10 works well for me
            open(myVideo)
            
            % change in tip deflection
            tip_0 = tip(1,:);
            tip_delta = tip(2:end, :) - tip(1:end-1, :);

            % Get the plot bounds for
            minZ = 1.2*min(min(tip));
            maxZ = 1.2*max(max(tip));
            
            % create the new figue
            %myVideoFigHandle = figure();
            myVideoFigHandle = figure('Renderer', 'opengl', 'Position', [10 10 900 800]);
            
            % Loop through all time steps
            for lv1 = 1:length(t)
            
                % Get current theta and deflection
                if lv1 ==1 % Nominal tip deflection
                    tip_i = tip_0;
                else % Update tip deflection to a new value with the delta
                    tip_i = tip_delta(lv1-1,:);
                end
            
                % Pause to let the figure catch up
                pause(1/100); %Pause and grab frame
            
                % Apply the new tip deflection 
                obj = obj.applyTipDeflection(tip_i(1), tip_i(2), tip_i(3), tip_i(4));
            
                %
                hold off
                obj.plotSailMeshOnFigure(z_scaling, myVideoFigHandle)
                % ax = gca;
                % ax.XAxis.Visible = 'off';
                % ax.YAxis.Visible = 'off';
                % ax.ZAxis.Visible = 'off';
                % gcf
                zlim([minZ, maxZ]);
                grid off
                pbaspect([1 1 0.4])  

                ax = gca;
                ax.FontSize = 16;
                fontSizeLabel = 26;
                xlabel('$b_1$ - meters', 'Interpreter','latex', 'FontSize',fontSizeLabel)
                ylabel('$b_2$ - meters', 'Interpreter','latex', 'FontSize',fontSizeLabel)
                zlabel('$b_3$ - meters', 'Interpreter','latex', 'FontSize',fontSizeLabel)
                set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
                pause(1/100); %Pause and grab frame
                %set(findall(myVideoFigHandle,'-property','FontSize'),'FontSize',26)
                %set(gcf,'units','points','position',[100,100,1000,800])
                %set(myVideoFigHandle,'renderer','Painters')

                hold off
            
                pause(1/100); %Pause and grab frame
                frame = getframe(myVideoFigHandle); %get frame
                writeVideo(myVideo, frame);
            
            end
            close(myVideo)
        end %FUNCTION

        function plotSailMeshOnFigure(obj, Zscale, figHandle)
            % Transparancy of edge
            edgeTransparancy = 0.1;
            
            % Call the designated figure to make it the current one
            figure(figHandle)

            % Plot quadrant 1
            s = trisurf(obj.sail1.T, obj.sail1.X, obj.sail1.Y, Zscale*obj.sail1.Z);
            s.EdgeAlpha = edgeTransparancy;
            % Plot quadrant 2
            hold on
            s = trisurf(obj.sail2.T, obj.sail2.X, obj.sail2.Y, Zscale*obj.sail2.Z);
            s.EdgeAlpha = edgeTransparancy;
            % Plot quadrant 3
            hold on
            s = trisurf(obj.sail3.T, obj.sail3.X, obj.sail3.Y, Zscale*obj.sail3.Z);
            s.EdgeAlpha = edgeTransparancy;
            % PLot quadrant 4
            hold on
            s = trisurf(obj.sail4.T, obj.sail4.X, obj.sail4.Y, Zscale*obj.sail4.Z);
            s.EdgeAlpha = edgeTransparancy;
            hold on  

            % Plot the booms as straight pipes
            z_BoomPlot = linspace(0,1,100);
            z_BoomPlot = z_BoomPlot.^2;
            plot3( linspace(obj.center(1),obj.tip1(1),100), linspace(obj.center(2), obj.tip1(2),100), z_BoomPlot*obj.tip1(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip2(1),100), linspace(obj.center(2), obj.tip2(2),100), z_BoomPlot*obj.tip2(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip3(1),100), linspace(obj.center(2), obj.tip3(2),100), z_BoomPlot*obj.tip3(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip4(1),100), linspace(obj.center(2), obj.tip4(2),100), z_BoomPlot*obj.tip4(3), 'k-', "LineWidth", 5 )
            
            
            % pbaspect([1 1 0.4])  
            % title(obj.name)
            % xlabel('$b_1$ - meters')
            % ylabel('$b_2$ - meters')
            % zlabel('$b_3$ - meters')
            % set(findall(gcf,'-property','FontSize'),'FontSize',26)
            % %set(findall(gcf,'-property','FontSize'),'FontSize',22)
        end % FUNCTION

        function plotSailMesh(obj, Zscale)
            % Transparancy of edge
            edgeTransparancy = 0.1;


            figure()
            s = trisurf(obj.sail1.T, obj.sail1.X, obj.sail1.Y, Zscale*obj.sail1.Z);
            s.EdgeAlpha = edgeTransparancy;
            hold on
            s = trisurf(obj.sail2.T, obj.sail2.X, obj.sail2.Y, Zscale*obj.sail2.Z);
            s.EdgeAlpha = edgeTransparancy;
            hold on
            s = trisurf(obj.sail3.T, obj.sail3.X, obj.sail3.Y, Zscale*obj.sail3.Z);
            s.EdgeAlpha = edgeTransparancy;
            hold on
            s = trisurf(obj.sail4.T, obj.sail4.X, obj.sail4.Y, Zscale*obj.sail4.Z);
            s.EdgeAlpha = edgeTransparancy;
            hold on

            z_BoomPlot = linspace(0,1,100);
            z_BoomPlot = z_BoomPlot.^2;
            plot3( linspace(obj.center(1),obj.tip1(1),100), linspace(obj.center(2), obj.tip1(2),100), z_BoomPlot*obj.tip1(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip2(1),100), linspace(obj.center(2), obj.tip2(2),100), z_BoomPlot*obj.tip2(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip3(1),100), linspace(obj.center(2), obj.tip3(2),100), z_BoomPlot*obj.tip3(3), 'k-', "LineWidth", 5 )
            plot3( linspace(obj.center(1),obj.tip4(1),100), linspace(obj.center(2), obj.tip4(2),100), z_BoomPlot*obj.tip4(3), 'k-', "LineWidth", 5 )
            pbaspect([1 1 0.4])  
            title(obj.name)

            xlabel('$b_1$ - meters')
            ylabel('$b_2$ - meters')
            zlabel('$b_3$ - meters')
            set(findall(gcf,'-property','FontSize'),'FontSize',26)
            %set(findall(gcf,'-property','FontSize'),'FontSize',22)
        end % FUNCTION

        function plotSRP_Torque(obj, SIA_span, CLOCK_span)
            figure()
            hold on
            sgtitle("Torque Total")
            lgd_str = cell(1, length(SIA_span));
            
            % Compute the torques
            torques = obj.computeTorqueAngles(SIA_span, CLOCK_span);
            

            for lv_sia = 1: length(SIA_span)
                subplot(3,1,1)
                hold on
                plot(CLOCK_span, torques(:,1,lv_sia) )
                lgd_str{lv_sia} = sprintf('%.0f',SIA_span(lv_sia));

                subplot(3,1,2)
                hold on
                plot(CLOCK_span, torques(:,2,lv_sia) )

                subplot(3,1,3)
                hold on
                plot(CLOCK_span, torques(:,3,lv_sia) )
            end
            
            subplot(3,1,1)
            if length(lgd_str)>1
                legend(lgd_str)
            end
            ylabel('$b_1$ Torque')
            xlim([0,360])
            title(obj.name)
            grid minor
            
            subplot(3,1,2)
            ylabel('$b_2$ Torque')
            grid minor
            xlim([0,360])
            
            subplot(3,1,3)
            ylabel('$b_3$ Torque')
            grid minor
            xlim([0,360])
            xlabel("Clock Angle - Degrees")
            set(findall(gcf,'-property','FontSize'),'FontSize',20)
        end % FUNCTION
        
    end % METHODS
    

    
% ----------------------------------------------------------------------- %
%                        DEPRACATED FUNCTIONS                             %
% ----------------------------------------------------------------------- %


    methods
        function [max1, max2, max3] = maxTorque(obj)
            % This should probably be a dynamic variable...
            % Find the maximum torque from the computed torques
    
            % Check if the torque is computed
            % if ~exist(obj.Torque_total_shape,'var')
            %     error('Torque not computed for this sail object.')
            % end

            % Get the torque value
            torqueValues = obj.Torque_total_shape;
            
            % Direction 1 ----------------------------------------------- %
            % obj.Torque_total_shape = zeros(3,length(CLOCK_span),length(SIA_span));
            max1_SIA = max(torqueValues(1, :, :), [], 2);

            % Overall max in the 1 direction
            max1 = max(max1_SIA);

            % Direction 2 ----------------------------------------------- %
            % obj.Torque_total_shape = zeros(3,length(CLOCK_span),length(SIA_span));
            max2_SIA = max(torqueValues(2, :, :), [], 2);

            % Overall max in the 2 direction
            max2 = max(max2_SIA);

            % Direction 3 ----------------------------------------------- %
            % obj.Torque_total_shape = zeros(3,length(CLOCK_span),length(SIA_span));
            max3_SIA = max(torqueValues(3, :, :), [], 2);

            % Overall max in the 3 direction
            max3 = max(max3_SIA);

        end

        function plotSRP_Force(obj)
            figure()
            hold on
            sgtitle("Force Total")
            lgd_str = cell(1, length(obj.SIA_span));
            
            for lv_sia = 1: length(obj.SIA_span)            
                subplot(3,1,1)
                hold on
                plot(obj.CLOCK_span, obj.Force_total_shape(1,:,lv_sia) )
                subplot(3,1,2)
                hold on
                plot(obj.CLOCK_span, obj.Force_total_shape(2,:,lv_sia) )
                subplot(3,1,3)
                hold on
                plot(obj.CLOCK_span, obj.Force_total_shape(3,:,lv_sia) )
            
                lgd_str{lv_sia} = sprintf('%.0f',obj.SIA_span(lv_sia));
            end
            
            subplot(3,1,1)
            legend(lgd_str)
            ylabel('$b_1$ Force')
            title(obj.name)
            xlim([0,360])
            
            grid minor
            subplot(3,1,2)
            ylabel('$b_2$ Force')
            grid minor
            xlim([0,360])
            
            subplot(3,1,3)
            ylabel('$b_3$ Force')
            grid minor
            xlim([0,360])
            xlabel("Clock Angle - Degrees")
            
            set(findall(gcf,'-property','FontSize'),'FontSize',20)
        end

        function torques = computeTorqueAngles(obj, SIA_span, CLOCK_span)
            % Compare torque and force generated vs flat plate
            Sun_dir_sun_frame = [0;0;1];

            torques = zeros(length(CLOCK_span), 3, length(SIA_span));
            
            % Inialize force and torque matrices
            obj.Torque_total_shape = zeros(3,length(CLOCK_span),length(SIA_span));
            
            for lv_sia = 1: length(SIA_span)
                for lv_clock = 1:length(CLOCK_span)
                    SIA = SIA_span(lv_sia);
                    clock_angle = CLOCK_span(lv_clock);
            
                    C_b_clock = [cosd(clock_angle) -sind(clock_angle) 0;
                                sind(clock_angle) cosd(clock_angle) 0;
                                 0 0 1 ];
                    C_clock_sun = [cosd(SIA) 0 sind(SIA);
                              0 1 0;
                              -sind(SIA) 0 cosd(SIA)];
            
                    C_a_sun = C_b_clock * C_clock_sun;
                    % assuming C_ab = eye(3) at t=0
                    Sun_dir_a = C_a_sun*Sun_dir_sun_frame;
            
                    % Form C_ab from time history
                    C_ab = eye(3);
            
                    % Transform sun into the body frame
                    Sun_dir_b = C_ab' * Sun_dir_a;
            
                    if ~isreal(obj.Torque_total_plate) || ~isreal(obj.Force_total_plate)
                        disp("Complex Torque or Force")
                    end

                    % Compute torques contributed from each sail quadrant
                    [~,Torque1] = MeshSailModel(obj.sail1.X, obj.sail1.Y, obj.sail1.Z, obj.sail1.T, Sun_dir_b, obj.Opt_prop);
                    [~,Torque2] = MeshSailModel(obj.sail2.X, obj.sail2.Y, obj.sail2.Z, obj.sail2.T, Sun_dir_b, obj.Opt_prop);
                    [~,Torque3] = MeshSailModel(obj.sail3.X, obj.sail3.Y, obj.sail3.Z, obj.sail3.T, Sun_dir_b, obj.Opt_prop);
                    [~,Torque4] = MeshSailModel(obj.sail4.X, obj.sail4.Y, obj.sail4.Z, obj.sail4.T, Sun_dir_b, obj.Opt_prop);

                    % Add torques to the results matrix
                    torques(lv_clock, :, lv_sia) = Torque1 + Torque2 + Torque3 + Torque4;
                end % END clock angle loop
            end % END SIA loop
        end % END Method function
  
    end %METHODS
end %CLASS












