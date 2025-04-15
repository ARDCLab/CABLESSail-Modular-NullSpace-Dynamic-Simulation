function Gamma_ia = compute_gamma(obj, t, coords, vels, q, v)
    %compute_gamma Summary of this method goes here
    %   Detailed explanation goes here
    %
    % ----------------------------------------------------------- %
    p_ia = coords{2};
    gamma_i = compute_gamma_vectorizedDCM(p_ia);
    Gamma_ia = blkdiag(eye(3), gamma_i, eye(length(coords{3})) );

    % Create base versions of these for converting DCM rates, euler angle
    % rate, etc to angular velocity. So you can either copy and paste that
    % function here, or just call it here and have those saved in the
    % utilities folder so the user can just use those.

end