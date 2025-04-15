function M = compute_M(obj, t, coords, vels, q, v)
    M = blkdiag(obj.const.mass_AMT*eye(3), obj.const.I_Bbb_b);
end