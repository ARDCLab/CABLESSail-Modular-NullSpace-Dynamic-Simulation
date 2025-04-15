function z = BasisGauvainTyler(x, y, L)
% This is the sail membrane basis function based on the 2023 ISSS paper
% from Andrew Gauvain and Daniel tyler from NASA Marshall Space Flight
% Center.
% A Solar Sail Shape Modeling Approach for Attitude Control Design and Analysis
%
% It is composed of a 2D sinusoid that matches the boundary condions of
% flat edges along the booms. The 2D centroid of the membrane is the high
% point of the basis function, which lies at (x, y) = (L/3, L/3), or 2/3
% along the line connecting the sail vertex to the center of the membrane
% edge. The basis function will always be <= 1, so the scaling that
% corresponds to the max z deflection in the paper will be handeled outside
% of this basis function in the sail mesh object method, similar to all
% other basis functions.
    r = sqrt(x^2+y^2);
    theta = atan2(y,x);
    z = sin(3*sqrt(2)*r*pi/(4*L)) * sin(2*theta);
end