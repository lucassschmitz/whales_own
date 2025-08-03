function L = globalLik_v2(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk)
% globalLik  Computes full-sample likelihood by multiplying captainLik
 
    C = max(c_id);
    C = 284; 
    L = 0; % initialize the log likelihood. 
    
    for c = 1:C
        mask  = (c_id == c);
        d_c    = d_vec(mask);
        Y_c    = Y_vec(mask);
        X_c    = X_mat(mask,:);
        Tau_c  = Tau_vec(mask);
        L_c    = captainLik(theta, d_c, Y_c, X_c, Tau_c, xk,wk);
        L      = L + log(L_c);
    end
end