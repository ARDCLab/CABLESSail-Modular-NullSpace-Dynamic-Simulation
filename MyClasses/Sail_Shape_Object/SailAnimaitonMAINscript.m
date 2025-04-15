% Animate sail shape with tip deflection from the scictech results

% Load the tip data from scitech resutls
load("SailVideoData.mat");
addpath('/Users/keeganbunker/Library/CloudStorage/GoogleDrive-bunke029@umn.edu/My Drive/Research-Work/Cable sail work/Sail shape generation/Sail_Shape_Object')

% Check the data --------------------------------------------------------
        % figure()
        % plot(time, b1_tip)
        % hold on
        % plot(time, b2_tip)
        % plot(time, b3_tip)
        % plot(time, b4_tip)


% Down sample the low freq data ------------------------------------------
index_select = 1;
b1Plot = b1_tip(1:index_select:end);
b2Plot = b2_tip(1:index_select:end);
b3Plot = b3_tip(1:index_select:end);
b4Plot = b4_tip(1:index_select:end);
timeplot = time(1:index_select:end);
        % figure()
        % plot(b1Plot)
        % hold on
        % plot(b2Plot)
        % plot(b3Plot)
        % plot(b4Plot)

% Create a sail ----------------------------------------------------------
% Sail Parameters
L_boom                      = 30;% Boom Length
n_nodes                     = 30;% Number of nodes
SailDeflectionAmplitude    = -0.5;% meters,
basisFunction               = @(x,y) BasisGauvainTyler(x, y, L_boom);

% Create sail object and apply tip deflection
sail_nominal = SailMesh(L_boom, n_nodes);
sail_nominal = sail_nominal.applyTipDeflection(0,0,0,0);
sail_nominal = sail_nominal.applyBasisFunction(1, basisFunction, SailDeflectionAmplitude);
sail_nominal = sail_nominal.applyBasisFunction(2, basisFunction, SailDeflectionAmplitude);
sail_nominal = sail_nominal.applyBasisFunction(3, basisFunction, SailDeflectionAmplitude);
sail_nominal = sail_nominal.applyBasisFunction(4, basisFunction, SailDeflectionAmplitude);
sail_nominal = sail_nominal.applyTipDeflection(0,0,0,0);

% Now animate! -----------------------------------------------------------
tipDeflections = [b1Plot', b2Plot', b3Plot', b4Plot'];
sail_nominal.AnimateSail(timeplot, tipDeflections, 1, pwd, 'Videotest')




