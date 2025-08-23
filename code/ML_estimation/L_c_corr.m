
function logL = L_c_corr(theta, a_c, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
% PROBLEM: this function does not use the log-sum-exp trick, causing underflowing issues. 
% 
% L_c_corr  Captain-level log-likelihood conditional on quality a_c
%   theta        full parameter vector (including product betas, alphas,
%                deltas, gamma0, gamma1, and covariance elements of ω)
%   a_c          scalar captain quality draw
%   d_cap        [J·Ncap×1] zero indicators for each product-voyage
%   Y_cap        [J·Ncap×1] observed outputs (Y>0 when d=1)
%   Xmat_cap     [J·Ncap×p] covariates for each product-voyage
%   tau_v_cap    [J·Ncap×1] scale parameters τ for each voyage (repeated J times)
%   xk, wk       1D Gauss–Hermite nodes & weights for a_c (not used here but
%                passed for consistency)
%   xk2,wk2,xk3,wk3  nodes & weights for multivariate GH in L_v_corr_v2
%
% Returns:
%   log_Lc       scalar log-likelihood of all voyages for this captain,
%                evaluated at a_c (to be integrated over a_c later)

    J     = 3;
    Ncap  = numel(d_c) / J;
    logL  = 0;

    logL = 0 ; % initialize log-likelihood 
    % Loop over each voyage for this captain
    for v = 1:Ncap
        idx    = (v-1)*J + (1:J);
        d_v    = d_c(idx);
        Y_v    = Y_c(idx);
        X_v    = x_c(idx(1), :);
        tau_v  = tau_c(idx(1));

        % Voyage-level log-likelihood with correlated ω (conditional on a_c)
        L_v = L_v_corr(theta, a_c, d_v, Y_v, X_v, tau_v, ...
                            xk, wk, xk2, wk2, xk3, wk3);

        logL = logL + log(L_v);
    end
end
