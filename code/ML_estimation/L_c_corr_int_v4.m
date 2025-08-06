
function log_Lc = L_c_corr_int_v4(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
    % difference with v3 is that uses the L_v_corr_v4 which introduces lambda parameters 
    J     = 3;
    Ncap  = numel(d_c) / J;
    Mq    = length(xk);
    
    ac_nodes = sqrt(2) * xk;
    log_w    = log(wk(:)) - 0.5*log(pi);
    log_tk   = zeros(Mq,1); % Initialize with zeros

    for qi = 1:Mq
        a_c   = ac_nodes(qi);
        logL_voyages  = 0; % Log-likelihood for all voyages of this captain for a given a_c

        for v = 1:Ncap
            idx   = (v-1)*J + (1:J);
            d_v   = d_c(idx);
            Y_v   = Y_c(idx);
            x_v   = x_c(idx(1), :);
            tau_v = tau_c(idx(1));
            
            logL_voyages = logL_voyages + L_v_corr_v4(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);
        end
        log_tk(qi) = log_w(qi) + logL_voyages;
    end

    % Log-sum-exp to compute final integrated log-likelihood
    M = max(log_tk);
    log_Lc = M + log(sum(exp(log_tk - M)));
end