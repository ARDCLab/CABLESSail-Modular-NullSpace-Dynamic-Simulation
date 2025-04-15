% Creates skey-symmetric 
%
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: August 15, 2023
%
%
%-------------------------------------------------------------------------%

function skew = crossr(r)

if size(r,1) ~= 3
    error("r must be a 3xn vector or matrix.")
end

n = size(r,2);
z = zeros(1,n);

skew = [z -r(3,:) r(2,:);
        r(3,:) z -r(1,:);
        -r(2,:) r(1,:) z];


% TODO: Add checks to make sure r is a vector. Throw warnings if a row
% vector is inputted.

