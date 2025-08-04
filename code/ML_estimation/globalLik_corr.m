function L = globalLik_corr(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk, xk2, wk2, xk3, wk3)
% globalLik_corr  Full-sample log-likelihood across all captains allowing
% for correlation. 
%   theta    parameter vector as above
%   d_vec,Y_vec,X_mat,Tau_vec  data stacks [J·Ntotal×1 or ×p]
%   c_id     [J·Ntotal×1] captain IDs (1...C)
%   xk,wk,...  GH nodes & weights

    C = max(c_id);
    logL_c = zeros(C,1);

    % ——— START parallel pool if possible ———
    try
        poolobj = gcp;
        if isempty(poolobj)
            parpool;
        end
    catch
        if exist('matlabpool','file') && matlabpool('size')==0
            matlabpool open;
        end
    end

    % Compute each captain's log-likelihood in parallel
    parfor c = 1:C
        mask   = (c_id == c);
        d_c    = d_vec(mask);
        Y_c    = Y_vec(mask);
        X_c    = X_mat(mask, :);
        Tau_c  = Tau_vec(mask);

        logL_c(c) = L_c_corr_int_v2(theta, d_c, Y_c, X_c, Tau_c, ...
                                    xk, wk, xk2, wk2, xk3, wk3);
    end

    % Sum across captains
    L = sum(logL_c);
end
