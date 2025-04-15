function ue = q2ue(x, q, psi)
% Compute the elastic coordinates, q, for all deflections, ue, at
% locations, x, with basis functions, psi.
%
% Inputs
%   x -     1D vector of x locations of length n
%   q -    n deflections
%   psi -   m basis functions
%
% Outputs
%   ue -     Deflections - n
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: May 31, 2023
%-------------------------------------------------------------------------%
% TODO:
%   -   Vectorize properly to compute multiple sets of ue at once with
%   proper matrix inputs for x and q

%-------------------------------------------------------------------------%

% Force vectors to be in the expected form
x = reshape(x, length(x),1);
q = reshape(q, length(q),1);

% compute deflections 
psiFT = matlabFunction(psi);
ue = psiFT(x) * q;
end