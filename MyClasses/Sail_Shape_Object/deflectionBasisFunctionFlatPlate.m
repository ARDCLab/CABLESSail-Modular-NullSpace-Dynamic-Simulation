function planeZ = deflectionBasisFunctionFlatPlate(tip1, tip2)
    center = [0 0 0]';
    % Sail normal vector
    planeNormal = cross( (tip2 - center), (tip1 - center));
    planeNormal = planeNormal/norm(planeNormal);
    a = planeNormal(1);
    b = planeNormal(2);
    c = planeNormal(3);
    
    % Apply plane deflection
    planeZ = @(x, y) -a/c*x - b/c*y;
end