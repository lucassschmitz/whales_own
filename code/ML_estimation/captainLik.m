
function logLc = captainLik(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk)
    % CAPTAIN-LEVEL likelihood using Gauss-Hermite with product-specific betas
    % and using logs to avoid underflowing. Returns the log captain
    % likelihood .
    J       = 3;
    Ncap    = numel(d_cap) / J;
    ac_nodes = sqrt(2) * xk;  % nodes for a_c ~ N(0,1)

    % Reshape stacked vectors
    d_mat    = reshape(d_cap, J, Ncap);
    Y_mat    = reshape(Y_cap, J, Ncap);
    tau_vec  = tau_v_cap(1:J:end)';             % 1×Ncap
    p        = size(Xmat_cap,2);
    X3d      = reshape(Xmat_cap, J, Ncap, p);   % J×Ncap×p

    % Unpack parameters
    beta_vec    = theta(1 : J*p);              
    beta_mat    = reshape(beta_vec, J, p);      % J×p matrix  
    idx         = J*p;                          

    alpha       = theta(idx+1    : idx+J);
    delta       = theta(idx+J+1  : idx+2*J);
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);
    sigma_omega = theta(idx+2*J+3: idx+3*J+2);

    Mq      = length(xk);
    logL_vals = zeros(Mq,1);             %stores logs to avoid underflow. 

    % Loop over Gauss-Hermite nodes
    for i = 1:Mq
        a_c = ac_nodes(i);

        % compute means for each product & voyage
        mu_part = zeros(J, Ncap);
        for kk = 1:p
            mu_part = mu_part + bsxfun(@times, X3d(:,:,kk), beta_mat(:,kk));
        end
        mu_z_mat = mu_part + delta * a_c;          % J×Ncap
        sd_mat   = sigma_omega * ones(1, Ncap);

        % per-product×voyage likelihoods
        u_raw   = bsxfun(@rdivide, log(Y_mat), alpha);
        u_val   = bsxfun(@minus, u_raw, log(tau_vec));
        p_val   = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
        f_z     = normpdf(u_val, mu_z_mat, sd_mat);
        jac     = 1 ./ (alpha * ones(1,Ncap) .* Y_mat);
        L_j_mat = p_val .* f_z .* jac;

        % zero-output cells replacement
        zero_idx = (d_mat == 0);
        if any(zero_idx(:))
            nodes     = bsxfun(@plus, reshape(mu_z_mat,1,J,Ncap), ...
                               reshape(sqrt(2)*sd_mat,1,J,Ncap).*reshape(xk,Mq,1,1));
            integrand = 1 ./ (1 + exp(gamma0 - gamma1 * nodes));
            weighted  = bsxfun(@times, reshape(wk,Mq,1,1), integrand);
            int_vals  = squeeze(sum(weighted,1)) / sqrt(pi);
            L_j_mat(zero_idx) = 1 - int_vals(zero_idx);
        end

        %L_j_mat = max( L_j_mat, realmin ); % floor negative values. 

        % accumulate log-likelihood
        logL_vals(i) = sum(log(L_j_mat(:)));     
    end
   
    % integrate via log-sum-exp
    log_w = log(wk(:)) - 0.5*log(pi);            
    log_tk = log_w + logL_vals;                 
    M      = max(log_tk);                      
    sum_exp = sum(exp(log_tk - M));            
    logLc  = M + log(sum_exp);                
    %L_c    = exp(logLc);                       
end

