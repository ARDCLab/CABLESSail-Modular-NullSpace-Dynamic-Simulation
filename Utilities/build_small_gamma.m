% Build small gamma matrix
function gamma_i = build_small_gamma(p_ia)
    p1 = p_ia(1:3);
    p2 = p_ia(4:6);
    p3 = p_ia(7:9);
    gamma_i = [zeros(3,1), -p3, p2;
        p3, zeros(3,1), -p1;
        -p2, p1, zeros(3,1)];
end