function obj = compute_energies(obj)
%compute_energies Computes the kinetic, potential, and total energy of the
%full system.
%
%   At every time step the coordinates and velocities are extracted and the
%   ODE matrices are compute to calculate the energies with the mass and
%   stiffness matrices and the reduced set of coordinates and velocities.

% ----------------------------------------------------------------------- %

    obj.PE = zeros(size(obj.t));
    obj.KE = zeros(2, length(obj.t));

    for lv1 = 1:length(obj.t)
        nu_hat = obj.velocity.history(:,lv1);
        q = obj.state.history(:,lv1);
        yi = [q; nu_hat];
        ti = obj.t(lv1);

        % Cell array of current states
        current_states = obj.state.decompose(q);

        % Cell array of current velocities
        current_velocities = obj.velocity.decompose(nu_hat);

        % Compute current null space matrix
        Upsilon     = obj.compute_upsilon(current_states, current_velocities);
        Upsilon_dot = obj.compute_upsilon_dot(current_states, current_velocities);

        % Compute pre-null-space states
        nu = Upsilon * nu_hat;
        
        % THIS WILL BREAK RIGHT HERE UNTIL I FIX UPSILON
        y_cell = obj.create_y_vectors(q, nu, nu_hat);

        % Current state matrices
        K       = obj.compute_K(ti, y_cell);
        M_bar   = obj.compute_M_bar(ti, y_cell);
        f_non   = obj.compute_f_non(ti, y_cell);
        f       = obj.compute_f(ti, y_cell);
        gamma   = obj.compute_Gamma(ti, y_cell);
        
        % Compute energies
        obj.KE(1,lv1) = 1/2 * nu'*M_bar*nu;
        obj.PE(1,lv1) = 1/2 * q'*K*q;

    end
    obj.TE = obj.KE(1,:) + obj.PE;
    obj.TE = abs(obj.TE -obj.TE(1,1));


    figure()

    subplot(3,1,1)
    plot(obj.t, obj.PE)
    grid on
    % ylabel('Potential Energy', 'FontSize', 8)
    ylabel('$B_S(t)$', 'FontSize', 8)
    subplot(3,1,2)
    plot(obj.t, obj.KE(1,:))
    grid on
    hold on
    % ylabel('Kinetic Energy', 'FontSize', 8)
    ylabel('$T_{Sa/a}(t)$', 'FontSize', 8)
    subplot(3,1,3)
    plot(obj.t, obj.TE)
    grid on
    % ylabel('Change in Total Energy', 'FontSize', 8)
    ylabel('$| \Sigma E(t) - \Sigma E(t_0) |$', 'FontSize', 8)
end