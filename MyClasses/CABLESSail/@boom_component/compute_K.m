function K = compute_K(obj, t, coords, vels, q, v)
    %compute_K Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ------------------------------------------------------------------- %

    % The stiffness is constant, so it is stored in the const structure.
    K = obj.const.K;
end