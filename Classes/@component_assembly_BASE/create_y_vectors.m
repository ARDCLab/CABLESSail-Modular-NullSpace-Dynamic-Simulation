function y_cell = create_y_vectors(obj, q, nu, nu_hat)
%METHOD1 Summary of this method goes here
%   Detailed explanation goes here

% ----------------------------------------------------------------------- %
    y_cell = cell(1, length(obj.components));
    n_q = 0;
    n_nu = 0;
    for lv1 = 1:length(obj.components)
        n_q_i = sum(obj.components{lv1}.state.dimensions);
        n_nu_i = sum(obj.components{lv1}.velocity.dimensions);
        y_i = [q(1+n_q:n_q+n_q_i, end); nu(1+n_nu:n_nu+n_nu_i, end)];
        n_q = n_q + n_q_i;
        n_nu = n_nu + n_nu_i;

        y_cell{lv1} = y_i;
    end
end