
function log_Lc = L_c_corr_int_v2(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
% captainLik_full  Combined captain-level log-likelihood
%   Integrates over captain quality a_c and sums voyage-level logs
%   theta    full parameter vector (β, α, δ, γ0, γ1, Σ_ω entries)
%   d_c      [J·Ncap×1] zero indicators
%   Y_c      [J·Ncap×1] outputs
%   x_c      [J·Ncap×p] covariates
%   tau_c    [J·Ncap×1] scale parameters
%   xk,wk    1D Gauss–Hermite nodes & weights for a_c
%   xk2,wk2, xk3,wk3  nodes & weights for multivariate GH inside voyage

    J     = 3;
    Ncap  = numel(d_c) / J;
    Mq    = length(xk);

    % Precompute transformed nodes and weight logs
    ac_nodes = sqrt(2) * xk;             % for a_c ~ N(0,1)
    log_w    = log(wk(:)) - 0.5*log(pi);
    log_tk   = nan(Mq,1);

    % Outer GH: integrate over a_c
    for qi = 1:Mq
        a_c   = ac_nodes(qi);
        logL  = 0;

        % Sum log-likelihood across voyages for this a_c
        for v = 1:Ncap
            idx   = (v-1)*J + (1:J);
            d_v   = d_c(idx);
            Y_v   = Y_c(idx);
            x_v   = x_c(idx(1), :);
            tau_v = tau_c(idx(1));

            % Voyage-level log-likelihood (multivariate GH)
            logL = logL + log(L_v_corr_v2(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3));
        end

        log_tk(qi) = log_w(qi) + logL;
    end

    % Log-sum-exp to compute final integrated log-likelihood
    M        = max(log_tk);
    log_Lc   = M + log(sum(exp(log_tk - M)));
end
