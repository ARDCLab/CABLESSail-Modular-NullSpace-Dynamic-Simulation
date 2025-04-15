% Build Large gamma matrix
function Gamma_ia = compute_gamma(obj, t, coords, vels, q, v)
    p_ia = coords{2};

    gamma_i = build_small_gamma(p_ia);
    Gamma_ia = blkdiag(eye(3),gamma_i);
end