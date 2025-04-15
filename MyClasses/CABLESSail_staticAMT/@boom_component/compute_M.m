function M = compute_M(obj, t, coords, vels, q, v)
    %compute_M Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ----------------------------------------------------------- %
    C_ia = p2DCM(coords{2});

    % Nominal Mass Matrix --------------------------------------- %
    M = ...
    [ obj.const.m * eye(3)  , -C_ia' * obj.const.c_B_ix , C_ia'*obj.const.H_B_i ;
    -obj.const.c_B_ix'*C_ia , obj.const.I_b_i           , -obj.const.G_B_i      ;
    obj.const.H_B_i'*C_ia   , -obj.const.G_B_i'         , obj.const.P_B_i ];
    % ----------------------------------------------------------- %

end