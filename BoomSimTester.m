% Boom tester

addpath(genpath([pwd, '/Utilities']))
addpath(genpath([pwd, '/Classes']))
set(0,'DefaultFigureWindowStyle','docked')

boom = boom_component();
tf = 10;
dtN = 100;
boom = boom.simulate_component(tf*10, dtN*10);

boom.coordinate.plot();