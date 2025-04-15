function q = ue2qe(x, ue, psi)
% Compute the elastic coordinates, q, for all deflections, ue, at
% locations, x, with basis functions, psi.
%
% Inputs
%   x -     1D vector of x locations of length n
%   ue -    n deflections
%   psi -   m basis functions
%
% Outputs
%   q -     Elastic coordinates - m
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: May 31, 2023
%-------------------------------------------------------------------------%
% TODO:
%   -   Vectorize properly to compute multiple sets of q at once with
%   proper matrix inputs for x and ue

%-------------------------------------------------------------------------%

ue = reshape(ue, length(ue),1);
x = reshape(x, length(x),1);
Psi = matlabFunction(psi);
A = Psi(x);
q = pinv(A)*ue; % Could use pinv or \ and it should be the same solution
end