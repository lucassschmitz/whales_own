
function log_Lc_int = L_c_int(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk, xk2, wk2, xk3, wk3)
% L_c_int  Captain-level log-likelihood integrated over a_c ~ N(0,1)
%   Integrates L_c_corr(theta, a_c) using 1D Gauss–Hermite quadrature
%
    % Number of GH nodes
    Mq = length(xk);
    logL_vals = zeros(Mq,1);

    % Transform nodes for standard normal
    ac_nodes = sqrt(2) * xk;

    % Evaluate conditional log-L at each node
    for i = 1:Mq
        logL_vals(i) = L_c_corr(theta, ac_nodes(i), ...
            d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk, xk2, wk2, xk3, wk3);
    end

    % Log-sum-exp to approximate the integral
    log_w  = log(wk(:)) - 0.5*log(pi);
    log_tk = log_w + logL_vals;
    M      = max(log_tk);
    log_Lc_int = M + log(sum(exp(log_tk - M)));
end
