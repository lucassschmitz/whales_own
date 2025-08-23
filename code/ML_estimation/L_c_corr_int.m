
function log_Lc = L_c_corr_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
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
            
            logL_voyages = logL_voyages + L_v_corr(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);
        end
        log_tk(qi) = log_w(qi) + logL_voyages;
    end

    % Log-sum-exp to compute final integrated log-likelihood
    M = max(log_tk);
    log_Lc = M + log(sum(exp(log_tk - M)));
end


%History 
% L_c_corr_int_v2 integrates directly over voyages whereas L_c_corr_int uses a L_c_corr. Hence L_c_corr multiplies voyages and then L_c_corr_int integrates over L_c_corr.
% L_c_corr_int_v3 improves on L_c_corr_int_v2 just changes the version of L_v_corr that uses from L_v_corr_v2 to L_v_corr_v3
% L_c_corr_int_v4 improves on L_c_corr_int_v3 just changes the version of L_v_corr that uses from L_v_corr_v3 to L_v_corr_v4

% function log_Lc_int = L_c_corr_int(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk, xk2, wk2, xk3, wk3)
% % L_c_corr is for a given captain skill, this function integrates over a_c.  
% %   Integrates L_c_corr(theta, a_c) using 1D Gauss–Hermite quadrature
% %
%     % Number of GH nodes
%     Mq = length(xk);
%     logL_vals = zeros(Mq,1);
% 
%     % Transform nodes for standard normal
%     ac_nodes = sqrt(2) * xk;
% 
%     % Evaluate conditional log-L at each node
%     for i = 1:Mq
%         logL_vals(i) = L_c_corr(theta, ac_nodes(i), ...
%             d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk, xk2, wk2, xk3, wk3);
%     end
% 
%     % Log-sum-exp to approximate the integral
%     log_w  = log(wk(:)) - 0.5*log(pi);
%     log_tk = log_w + logL_vals;
%     M      = max(log_tk);
%     log_Lc_int = M + log(sum(exp(log_tk - M)));
% end

% function log_Lc = L_c_corr_int_v2(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
% % captainLik_full  Combined captain-level log-likelihood
% %   Integrates over captain quality a_c and sums voyage-level logs
% %   theta    full parameter vector (β, α, δ, γ0, γ1, Σ_ω entries)
% %   d_c      [J·Ncap×1] zero indicators
% %   Y_c      [J·Ncap×1] outputs
% %   x_c      [J·Ncap×p] covariates
% %   tau_c    [J·Ncap×1] scale parameters
% %   xk,wk    1D Gauss–Hermite nodes & weights for a_c
% %   xk2,wk2, xk3,wk3  nodes & weights for multivariate GH inside voyage
% 
%     J     = 3;
%     Ncap  = numel(d_c) / J;
%     Mq    = length(xk);
% 
%     % Precompute transformed nodes and weight logs
%     ac_nodes = sqrt(2) * xk;             % for a_c ~ N(0,1)
%     log_w    = log(wk(:)) - 0.5*log(pi);
%     log_tk   = nan(Mq,1);
% 
%     % Outer GH: integrate over a_c
%     for qi = 1:Mq
%         a_c   = ac_nodes(qi);
%         logL  = 0;
% 
%         % Sum log-likelihood across voyages for this a_c
%         for v = 1:Ncap
%             idx   = (v-1)*J + (1:J);
%             d_v   = d_c(idx);
%             Y_v   = Y_c(idx);
%             x_v   = x_c(idx(1), :);
%             tau_v = tau_c(idx(1));
% 
%             % Voyage-level log-likelihood (multivariate GH)
%             logL = logL + log(L_v_corr_v2(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3));
%         end
% 
%         log_tk(qi) = log_w(qi) + logL;
%     end
% 
%     % Log-sum-exp to compute final integrated log-likelihood
%     M        = max(log_tk);
%     log_Lc   = M + log(sum(exp(log_tk - M)));
% end

% function log_Lc = L_c_corr_int_v3(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)
%     % --- Integrates over captain quality a_c and sums voyage-level LOGS ---
%     J     = 3;
%     Ncap  = numel(d_c) / J;
%     Mq    = length(xk);
% 
%     ac_nodes = sqrt(2) * xk;
%     log_w    = log(wk(:)) - 0.5*log(pi);
%     log_tk   = zeros(Mq,1); % Initialize with zeros
% 
%     for qi = 1:Mq
%         a_c   = ac_nodes(qi);
%         logL_voyages  = 0; % Log-likelihood for all voyages of this captain for a given a_c
% 
%         for v = 1:Ncap
%             idx   = (v-1)*J + (1:J);
%             d_v   = d_c(idx);
%             Y_v   = Y_c(idx);
%             x_v   = x_c(idx(1), :);
%             tau_v = tau_c(idx(1));
% 
%             logL_voyages = logL_voyages + L_v_corr_v3(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);
%         end
%         log_tk(qi) = log_w(qi) + logL_voyages;
%     end
% 
%     % Log-sum-exp to compute final integrated log-likelihood
%     M = max(log_tk);
%     log_Lc = M + log(sum(exp(log_tk - M)));
% end
 