function Mdot = compute_Mdot(obj, t, coords, vels, q, v)
    %compute_Mdot Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ----------------------------------------------------------- %
    [~, ~, C_ia_dot] = decompose_attitude(coords{2}, vels{2});
    
    Mdot = ...
    [ zeros(3,3)                , -C_ia_dot' * obj.const.c_B_ix , C_ia_dot'*obj.const.H_B_i ;
    -obj.const.c_B_ix'*C_ia_dot , zeros(3,3)                    , zeros(3,6);
    obj.const.H_B_i'*C_ia_dot   , zeros(6,3)                    , zeros(6,6)];
    
end