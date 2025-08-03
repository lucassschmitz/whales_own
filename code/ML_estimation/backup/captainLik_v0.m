
 
function L_c = captainLik(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk)
    % Vectorized captain-level likelihood using Gauss-Hermite
    J       = 3;
    Ncap    = numel(d_cap) / J;
    ac_nodes = sqrt(2) * xk;  % nodes for a_c ~ N(0,1)

    % Reshape stacked vectors into matrices: J rows × Ncap columns
    d_mat    = reshape(d_cap, J, Ncap);
    Y_mat    = reshape(Y_cap, J, Ncap);
    tau_vec  = tau_v_cap(1:J:end)';           % 1×Ncap row vector
    K        = size(Xmat_cap,2);
    X3d      = reshape(Xmat_cap, J, Ncap, K); % J×Ncap×K

    % Unpack parameters
    beta        = theta(1:2);
    alpha       = theta(3:5);
    delta       = theta(6:8);
    gamma0      = theta(9);
    gamma1      = theta(10);
    sigma_omega = theta(11:13);

    % Preallocate likelihoods per quadrature node
    Mq      = length(xk);
    L_vals  = zeros(Mq, 1);

    % Loop over Gauss-Hermite nodes
    for i = 1:Mq
        a_c = ac_nodes(i);

        % Linear predictor for all voyages at once
        mu_part = zeros(J, Ncap);
        for kk = 1:K
            mu_part = mu_part + X3d(:,:,kk) * beta(kk);
        end
        mu_z_mat = mu_part + delta * a_c;            % J×Ncap
        sd_mat   = sigma_omega * ones(1, Ncap);       % J×Ncap

        % Positive-output likelihood
        u_raw = bsxfun(@rdivide, log(Y_mat), alpha);             % J×Ncap
        u_val = bsxfun(@minus, u_raw, log(tau_vec));             % J×Ncap (tau_vec is 1×Ncap)
        p_val = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
        f_z   = normpdf(u_val, mu_z_mat, sd_mat);
        jac   = 1 ./ (alpha * ones(1, Ncap) .* Y_mat);
        L_j_mat = p_val .* f_z .* jac;

        % Zero-output adjustment via GH
        zero_idx = (d_mat == 0);
        if any(zero_idx(:))
            nodes     = bsxfun(@plus, reshape(mu_z_mat, 1, J, Ncap), ...
                               reshape(sqrt(2) * sd_mat, 1, J, Ncap) .* reshape(xk, Mq, 1, 1));
            integrand = 1 ./ (1 + exp(gamma0 - gamma1 * nodes));
            weighted  = bsxfun(@times, reshape(wk, Mq, 1, 1), integrand);
            int_vals  = squeeze(sum(weighted, 1)) / sqrt(pi);        % J×Ncap
            L_j_mat(zero_idx) = 1 - int_vals(zero_idx);
        end

        % Product over products and voyages
        L_vals(i) = prod(L_j_mat(:));
    end

    % Aggregate over nodes
    L_c = sum(wk .* L_vals) / sqrt(pi);
end


 