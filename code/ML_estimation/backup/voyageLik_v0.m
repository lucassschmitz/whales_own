function L_v = voyageLik(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk )
    J = 3;

    % Unpack parameters
    beta        = theta(1:2);
    alpha       = theta(3:5);
    delta       = theta(6:8);
    gamma0      = theta(9);
    gamma1      = theta(10);
    sigma_omega = theta(11:13);

    % Precompute: mu_z and sd_z for all j
    mu_z = x_v' * beta + delta * a_c;  % 3x1
    sd_z = sigma_omega;

    % Precompute common terms
    u_val  = log(Y_v) ./ alpha - log(tau_v);
    p_val  = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
    f_z    = normpdf(u_val, mu_z, sd_z);
    jac    = 1 ./ (alpha .* Y_v);

    % Vectorized likelihoods for d_v == 1 and d_v == 0
    L_j = zeros(J,1);

    pos_idx = d_v == 1;
    L_j(pos_idx) = p_val(pos_idx) .* f_z(pos_idx) .* jac(pos_idx);

    % Use Gauss-Hermite quadrature for the zero-output case
    for j = find(~pos_idx)'
        % Change of variables: u = sqrt(2)*sd*xk + mu
        nodes = sqrt(2) * sd_z(j) * xk + mu_z(j);
        weights = wk .* (1 ./ (1 + exp(gamma0 - gamma1 * nodes)));
        L_j(j) = 1 - sum(weights) / sqrt(pi);
    end

    % Final voyage likelihood is product across j
    L_v = prod(L_j);
end

