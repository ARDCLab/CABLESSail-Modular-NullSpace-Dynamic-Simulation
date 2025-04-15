function KE = KE_velocity(q_state, nu, boom)

syms x
Psi = boom.PSI;

% extract vectors
w = nu(4:6,1);
q_eps = q_state(13:end, 1);
q_eps_dot = nu(7:12, 1);
r_ca_dota_a = nu(1:3,1);
p_ia = q_state(4:12,1);
[C_ac, ~] = build_C_ab(p_ia, false);

% Compute positions and velocities
u_dot = Psi*q_eps_dot;
r_dm_c_dotc_c = [0;
               u_dot];

r_dm_c_c = [x
            Psi*q_eps];

r_dm_c_dota_c = crossr(w)*r_dm_c_c + r_dm_c_dotc_c;

r_ca_dota_c = C_ac' * r_ca_dota_a;

r_dma_dota_c = r_ca_dota_c + r_dm_c_dota_c;

e_fn = r_dma_dota_c'*r_dma_dota_c;


KE = 1/2*int(e_fn, 0, boom.L);
