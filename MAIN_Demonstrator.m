% Main file to run the open loop simulation of CABLESSail with the new
% somulation structure

%% Add folders to workspace --------------------------------------------- %
addpath(genpath([pwd, '/Utilities']))
addpath(genpath([pwd, '/Classes']))
set(0,'DefaultFigureWindowStyle','docked')

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Hub Component tester

% Create components
hub1 = hub_component();

% Simulate
tf = 10;
dtN = 100;
hub1 = hub1.simulate_component(tf, dtN);

%% Plot results
% Plot position
hub1.coordinate.states{1}.plot()

% Plot velocity
hub1.velocity.states{1}.plot()

% Plot attitude
hub1.coordinate.states{2}.plot()

% Plot angular velocity
hub1.velocity.states{2}.plot()

%% Boom component tester

% Create boom component
boom1 = boom_component();
 boom1.coordinate.states{3}.y0 = [1 0 0 0 0 0]';

% Simulate component
tf = 10;
dtN = 100;
boom1 = boom1.simulate_component(tf, dtN);

% Plot
% % Plot position
% boom1.coordinate.states{1}.plot()
% 
% % Plot velocity
% boom1.velocity.states{1}.plot()
% 
% % Plot attitude
% boom1.coordinate.states{2}.plot()
% 
% % Plot angular velocity
% boom1.velocity.states{2}.plot()

% Plot elastic coordinates
boom1.coordinate.states{3}.plot()











