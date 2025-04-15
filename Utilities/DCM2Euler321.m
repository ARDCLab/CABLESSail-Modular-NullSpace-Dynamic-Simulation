function [phi,theta,psi] = DCM2Euler321(Cba)
%DCM2EULER321  DCM to 3-2-1 Euler sequence angle conversion.
% [phi,theta,psi] = DCM2EULER321(Cba) solves for Euler angles associated 
% with 3-2-1 rotation sequence based on the DCM using procedure on p. 23 
% of de Ruiter (2013) 
%
% INPUT PARAMETERS:
% Cba = 3x3 DCM
%
% OUTPUT PARAMETERS:
% phi,theta,psi = Euler angles associated with 3-2-1 rotation sequence 
%                 calculated from Cba
%
% Ryan Caverly
% Updated February 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Euler angles
phi = atan2(Cba(2,3),Cba(3,3));
theta = -asin(Cba(1,3));
psi = atan2(Cba(1,2),Cba(1,1));
