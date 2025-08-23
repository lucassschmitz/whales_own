function L = globalLik(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk)
% globalLik  Computes full-sample likelihood by multiplying captainLik
%this version uses parfor to be more efficient by using parallel
%computation 
    C = max(c_id);  
    
    % preallocate one log‐likelihood per captain
    logL_c = zeros(C,1);

    % ——— START a pool if possible ———
    % Try the new gcp/no-create API; if that fails, fall back to matlabpool
    try
        % new interface (R2013b+)
        poolobj = gcp;         % no args
        if isempty(poolobj)
            parpool;          % launch default
        end
    catch
        % older interface
        if exist('matlabpool','file') && matlabpool('size')==0
            matlabpool open;
        end
    end
   
    parfor c = 1:C
        mask   = (c_id == c);
        d_c    = d_vec(mask);
        Y_c    = Y_vec(mask);
        X_c    = X_mat(mask,:);
        Tau_c  = Tau_vec(mask);
        % returns log L_c
        logL_c(c) = captainLik(theta, d_c, Y_c, X_c, Tau_c, xk, wk);
    end
    
    % sum up all the per‐captain logs
    L = sum(logL_c);
end
