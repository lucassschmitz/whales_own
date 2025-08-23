
function logL_v  = L_v_corr(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3)
    % the differ with v3 is that includes lambda. 
    
    
    J = numel(d_v);
    p = numel(x_v);
    
    %—— Unpack parameters
    beta_vec    = theta(1:J*p);
    beta_mat    = reshape(beta_vec, J, p);        % J×p
    idx         = J*p;
    alpha       = theta(idx+1    : idx+J);        % J×1
    delta       = theta(idx+J+1  : idx+2*J);      % J×1
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);
    Omega_vec   = theta(idx+2*J+3 : idx+2*J+2+J^2);
    Sigma_omega = reshape(Omega_vec, J, J);       % J×J covariance
    lambda = theta(end-J+1: end); 


    %—— Compute latent log-shocks
    w_hat = (log(Y_v) - log(lambda)) ./ alpha - log(tau_v);  % J×1
    %—— Partition indices
    pos_idx  = find(d_v==1);
    zero_idx = find(d_v==0);

    mu      = beta_mat * x_v + delta * a_c;       % J×1
    Sigma_pp= Sigma_omega(pos_idx,pos_idx);       % |J+|×|J+|
    mu_pos  = mu(pos_idx);                        % |J+|×1
    
    %—— Part 1: positive outputs (joint density)
    if isempty(pos_idx)
            %  no positive outputs ⇒ L_pos = 1
            log_L_pos = 0;   
    else
        
        % product of selection probabilities
        p_prod_log = 0;
        for ii = 1:numel(pos_idx)
            j = pos_idx(ii);
            log_p_j = -log(1 + exp(gamma0 - gamma1 * w_hat(j)));
            p_prod_log = p_prod_log + log_p_j; 
        end
        
        % joint density of w_hat over positives
        f_pos_joint = mvnpdf(w_hat(pos_idx)', mu_pos', Sigma_pp);  % scalar        
        jac_pos = prod( 1 ./ ( alpha(pos_idx) .* Y_v(pos_idx) ) );  
        log_L_pos = p_prod_log + log(f_pos_joint) + log(jac_pos);                  
    end

    %—— Part 2: zero outputs (one joint integral)
    % compute conditional mean and covariance
    mu_zero   = mu(zero_idx);
    Sigma_pp  = Sigma_omega(pos_idx,pos_idx);
    Sigma_zp  = Sigma_omega(zero_idx,pos_idx);
    Sigma_zz  = Sigma_omega(zero_idx,zero_idx);
    mu_cond   = mu_zero + Sigma_zp * (Sigma_pp \ (w_hat(pos_idx)-mu_pos));
    Sigma_cond= Sigma_zz - Sigma_zp * (Sigma_pp \ Sigma_omega(pos_idx,zero_idx));

    d0 = numel(zero_idx);
    if d0==0
        log_L_zero = 0; 
    elseif d0==1
        w_nodes = sqrt(2*Sigma_cond) * xk + mu_cond;
        f_nodes = exp(gamma0 - gamma1 .* w_nodes) ./ (1 + exp(gamma0 - gamma1 .* w_nodes));        
        L_zero = sum( wk .* f_nodes ) / sqrt(pi);
        log_L_zero = log(max(L_zero, realmin)); % avoid underflow

    elseif d0==2
        L = chol(Sigma_cond, 'lower');
        Z = sqrt(2) * (xk2 * L') + mu_cond'; 
        w_mv = wk2 / (pi);   
        
        logit_vals = 1 ./ (1 + exp( -(gamma0 - Z *gamma1) ));
        gvals      = prod(logit_vals, 2);          
        
        L_zero = sum(w_mv.* gvals); 
        log_L_zero =  log(max(L_zero, realmin)); % avoid underflow

    else % for d0>=3
        L = chol(Sigma_cond, 'lower');
        Z = sqrt(2) * (xk3 * L') + mu_cond'; 
        w_mv = wk3 / (pi^(3/2));   

        logit_vals = 1 ./ (1 + exp( -(gamma0 - Z * gamma1) ));
        gvals      = prod(logit_vals, 2);          

        L_zero = sum(w_mv.* gvals); 
        log_L_zero = log(max(L_zero, realmin));  % avoid underflow
    end

    logL_v = log_L_pos + log_L_zero; 
end


%History 
% L_v_corr_v2 produces same output as L_v_corr but using quadrature to make it faster. 
% L_v_corr_v3 takes L_v_corr_v2 but returnes the log-lik instead of lik, uses sum-log-exp trick. the log and the trick are to avoid underflow. 
% L_v_corr_v4 updates L_v_corr_v3 to allow lambda in the model, which multiplies production to adjust for the different units. 


% function L_v = L_v_corr(theta, a_c, d_v, Y_v, x_v, tau_v)
%     % implements voyage likelihood given with correlation  
%     J = numel(d_v);
%     p = numel(x_v);
% 
%     %—— Unpack parameters
%     beta_vec    = theta(1:J*p);
%     beta_mat    = reshape(beta_vec, J, p);        % J×p
%     idx         = J*p;
%     alpha       = theta(idx+1    : idx+J);        % J×1
%     delta       = theta(idx+J+1  : idx+2*J);      % J×1
%     gamma0      = theta(idx+2*J+1);
%     gamma1      = theta(idx+2*J+2);
%     Omega_vec   = theta(idx+2*J+3 : idx+2*J+2+J^2);
%     Sigma_omega = reshape(Omega_vec, J, J);       % J×J covariance
% 
%     %—— Compute latent log-shocks
%     w_hat = log(Y_v)./alpha - log(tau_v);       % J×1
% 
%     %—— Partition indices
%     pos_idx  = find(d_v==1);
%     zero_idx = find(d_v==0);
% 
% 
%     mu      = beta_mat * x_v + delta * a_c;       % J×1
%     Sigma_pp= Sigma_omega(pos_idx,pos_idx);       % |J+|×|J+|
%     mu_pos  = mu(pos_idx);                        % |J+|×1
% 
%     %—— Part 1: positive outputs (joint density)
%     if isempty(pos_idx)
%             % Change: no positive outputs ⇒ L_pos = 1
%             L_pos = 1;   
%     else
% 
%         % product of selection probabilities
%         p_prod = 1;
%         for ii = 1:numel(pos_idx)
%             j = pos_idx(ii);
%             p_j = 1/(1+exp(gamma0 - gamma1 * w_hat(j)));
%             p_prod = p_prod * p_j;
%         end
% 
%         % joint density of w_hat over positives
%         f_pos_joint = mvnpdf(w_hat(pos_idx)', mu_pos', Sigma_pp);  % scalar
% 
% 
%         % multiply by Jacobian ∏(1/(α_j·y_vj)) per eq. 18
%         jac_pos = prod( 1 ./ ( alpha(pos_idx) .* Y_v(pos_idx) ) );  
%         L_pos   = p_prod * f_pos_joint * jac_pos;                   
%     end
% 
%     %—— Part 2: zero outputs (one joint integral)
%     % compute conditional mean and covariance
%     mu_zero   = mu(zero_idx);
%     Sigma_pp  = Sigma_omega(pos_idx,pos_idx);
%     Sigma_zp  = Sigma_omega(zero_idx,pos_idx);
%     Sigma_zz  = Sigma_omega(zero_idx,zero_idx);
%     mu_cond   = mu_zero + Sigma_zp * (Sigma_pp \ (w_hat(pos_idx)-mu_pos));
%     Sigma_cond= Sigma_zz - Sigma_zp * (Sigma_pp \ Sigma_omega(pos_idx,zero_idx));
% 
%     d0 = numel(zero_idx);
%     if d0==0
%         L_zero = 1;
%     elseif d0==1
%         integrand = @(w) (1 ./ (1 + exp(-(gamma0 - gamma1.*w)))) .* normpdf(w, mu_cond, Sigma_cond);
%         L_zero = integral(integrand, -Inf, Inf);
% 
%     elseif d0==2
%         integrand2 = @(w1,w2) prod(exp(gamma0 - gamma1*[w1;w2])./(1+exp(gamma0 - gamma1*[w1;w2]))) .* ...
%                   mvnpdf([w1,w2], mu_cond', Sigma_cond);
%         L_zero = integral2(integrand2, -Inf, Inf, -Inf, Inf);
%     else % for d0>=3
% 
%         %integrand3 = @(w1,w2,w3) prod(exp(gamma0 - gamma1*[w1;w2;w3])./(1+exp(gamma0 - gamma1*[w1;w2;w3]))) .* ...
%         %              mvnpdf([w1,w2,w3], mu_cond', Sigma_cond);
% 
%         % integrand3 = @(w1,w2,w3) arrayfun( ...
%         %   @(z1,z2,z3) ...
%         %     prod( 1 ./ (1 + exp-((gamma0 - gamma1 .* [z1; z2; z3]))) ) ...
%         %   * mvnpdf([z1, z2, z3], mu_cond', Sigma_cond), ... 
%         %   w1, w2, w3);    % Change: wrap in arrayfun to force scalar mvnpdf calls
% 
%         integrand3 = @(w1,w2,w3) arrayfun( ...
%       @(z1,z2,z3) ...
%         prod( 1 ./ (1 + exp( -(gamma0 - gamma1 .* [z1; z2; z3]) ) ) ) ...
%       .* mvnpdf([z1, z2, z3], mu_cond', Sigma_cond), ... 
%       w1, w2, w3);
% 
%         L_zero = integral3( integrand3, -Inf, Inf, -Inf, Inf, -Inf, Inf);
%         %L_zero = integral3( integrand3, -Inf, Inf, -Inf, Inf, -Inf, Inf, 'ArrayValued', true );
%         %L_zero = integral(@(w1) integral(@(w2) integral(@(w3) integrand3(w1,w2,w3), -Inf, Inf), -Inf, Inf), -Inf, Inf);
%     end
% 
%     L_v = L_pos * L_zero;
% end

% function L_v  = L_v_corr_v2(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3)
%     % implements L_v_corr (produces same result) but trying to optimize the integration,
%     % spetially when the voyage produces 0 output the 3D integral is hard to compute. 
%     % implements voyage likelihood (NOT log-likelihood) given with correlation  
%     J = numel(d_v);
%     p = numel(x_v);
% 
%     %—— Unpack parameters
%     beta_vec    = theta(1:J*p);
%     beta_mat    = reshape(beta_vec, J, p);        % J×p
%     idx         = J*p;
%     alpha       = theta(idx+1    : idx+J);        % J×1
%     delta       = theta(idx+J+1  : idx+2*J);      % J×1
%     gamma0      = theta(idx+2*J+1);
%     gamma1      = theta(idx+2*J+2);
%     Omega_vec   = theta(idx+2*J+3 : idx+2*J+2+J^2);
%     Sigma_omega = reshape(Omega_vec, J, J);       % J×J covariance
% 
%     %—— Compute latent log-shocks
%     w_hat = log(Y_v)./alpha - log(tau_v);       % J×1
% 
%     %—— Partition indices
%     pos_idx  = find(d_v==1);
%     zero_idx = find(d_v==0);
% 
%     mu      = beta_mat * x_v + delta * a_c;       % J×1
%     Sigma_pp= Sigma_omega(pos_idx,pos_idx);       % |J+|×|J+|
%     mu_pos  = mu(pos_idx);                        % |J+|×1
% 
%     %—— Part 1: positive outputs (joint density)
%     if isempty(pos_idx)
%             %  no positive outputs ⇒ L_pos = 1
%             L_pos = 1;   
%     else
% 
%         % product of selection probabilities
%         p_prod = 1;
%         for ii = 1:numel(pos_idx)
%             j = pos_idx(ii);
%             p_j = 1/(1+exp(gamma0 - gamma1 * w_hat(j)));
%             p_prod = p_prod * p_j;
%         end
% 
%         % joint density of w_hat over positives
%         f_pos_joint = mvnpdf(w_hat(pos_idx)', mu_pos', Sigma_pp);  % scalar        
%         jac_pos = prod( 1 ./ ( alpha(pos_idx) .* Y_v(pos_idx) ) );  
%         L_pos   = p_prod * f_pos_joint * jac_pos;                   
%     end
% 
%     %—— Part 2: zero outputs (one joint integral)
%     % compute conditional mean and covariance
%     mu_zero   = mu(zero_idx);
%     Sigma_pp  = Sigma_omega(pos_idx,pos_idx);
%     Sigma_zp  = Sigma_omega(zero_idx,pos_idx);
%     Sigma_zz  = Sigma_omega(zero_idx,zero_idx);
%     mu_cond   = mu_zero + Sigma_zp * (Sigma_pp \ (w_hat(pos_idx)-mu_pos));
%     Sigma_cond= Sigma_zz - Sigma_zp * (Sigma_pp \ Sigma_omega(pos_idx,zero_idx));
% 
%     d0 = numel(zero_idx);
%     if d0==0
%         L_zero = 1;
%     elseif d0==1
%         w_nodes = sqrt(2*Sigma_cond) * xk + mu_cond;
%         f_nodes = exp(gamma0 - gamma1 .* w_nodes) ./ (1 + exp(gamma0 - gamma1 .* w_nodes));        
%         L_zero = sum( wk .* f_nodes ) / sqrt(pi);
% 
%     elseif d0==2
%         L = chol(Sigma_cond, 'lower');
%         Z = sqrt(2) * (xk2 * L') + mu_cond'; 
%         w_mv = wk2 / (pi);   
% 
%         logit_vals = 1 ./ (1 + exp( -(gamma0 - Z *gamma1) ));
%         gvals      = prod(logit_vals, 2);          
% 
%         L_zero = sum(w_mv.* gvals); 
% 
%     else % for d0>=3
%         L = chol(Sigma_cond, 'lower');
%         Z = sqrt(2) * (xk3 * L') + mu_cond'; 
%         w_mv = wk3 / (pi^(3/2));   
% 
%         logit_vals = 1 ./ (1 + exp( -(gamma0 - Z * gamma1) ));
%         gvals      = prod(logit_vals, 2);          
% 
%         L_zero = sum(w_mv.* gvals); 
%     end
% 
%     L_v = L_pos * L_zero;
% end

% function logL_v  = L_v_corr_v3(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3)
% 
%     J = numel(d_v);
%     p = numel(x_v);
% 
%     %—— Unpack parameters
%     beta_vec    = theta(1:J*p);
%     beta_mat    = reshape(beta_vec, J, p);        % J×p
%     idx         = J*p;
%     alpha       = theta(idx+1    : idx+J);        % J×1
%     delta       = theta(idx+J+1  : idx+2*J);      % J×1
%     gamma0      = theta(idx+2*J+1);
%     gamma1      = theta(idx+2*J+2);
%     Omega_vec   = theta(idx+2*J+3 : idx+2*J+2+J^2);
%     Sigma_omega = reshape(Omega_vec, J, J);       % J×J covariance
% 
%     %—— Compute latent log-shocks
%     w_hat = log(Y_v)./alpha - log(tau_v);       % J×1
% 
%     %—— Partition indices
%     pos_idx  = find(d_v==1);
%     zero_idx = find(d_v==0);
% 
%     mu      = beta_mat * x_v + delta * a_c;       % J×1
%     Sigma_pp= Sigma_omega(pos_idx,pos_idx);       % |J+|×|J+|
%     mu_pos  = mu(pos_idx);                        % |J+|×1
% 
%     %—— Part 1: positive outputs (joint density)
%     if isempty(pos_idx)
%             %  no positive outputs ⇒ L_pos = 1
%             log_L_pos = 0;   
%     else
% 
%         % product of selection probabilities
%         p_prod_log = 0;
%         for ii = 1:numel(pos_idx)
%             j = pos_idx(ii);
%             log_p_j = -log(1 + exp(gamma0 - gamma1 * w_hat(j)));
%             p_prod_log = p_prod_log + log_p_j; 
%         end
% 
%         % joint density of w_hat over positives
%         f_pos_joint = mvnpdf(w_hat(pos_idx)', mu_pos', Sigma_pp);  % scalar        
%         jac_pos = prod( 1 ./ ( alpha(pos_idx) .* Y_v(pos_idx) ) );  
%         log_L_pos = p_prod_log + log(f_pos_joint) + log(jac_pos);                  
%     end
% 
%     %—— Part 2: zero outputs (one joint integral)
%     % compute conditional mean and covariance
%     mu_zero   = mu(zero_idx);
%     Sigma_pp  = Sigma_omega(pos_idx,pos_idx);
%     Sigma_zp  = Sigma_omega(zero_idx,pos_idx);
%     Sigma_zz  = Sigma_omega(zero_idx,zero_idx);
%     mu_cond   = mu_zero + Sigma_zp * (Sigma_pp \ (w_hat(pos_idx)-mu_pos));
%     Sigma_cond= Sigma_zz - Sigma_zp * (Sigma_pp \ Sigma_omega(pos_idx,zero_idx));
% 
%     d0 = numel(zero_idx);
%     if d0==0
%         log_L_zero = 0; 
%     elseif d0==1
%         w_nodes = sqrt(2*Sigma_cond) * xk + mu_cond;
%         f_nodes = exp(gamma0 - gamma1 .* w_nodes) ./ (1 + exp(gamma0 - gamma1 .* w_nodes));        
%         L_zero = sum( wk .* f_nodes ) / sqrt(pi);
%         log_L_zero = log(max(L_zero, realmin)); % avoid underflow
% 
%     elseif d0==2
%         L = chol(Sigma_cond, 'lower');
%         Z = sqrt(2) * (xk2 * L') + mu_cond'; 
%         w_mv = wk2 / (pi);   
% 
%         logit_vals = 1 ./ (1 + exp( -(gamma0 - Z *gamma1) ));
%         gvals      = prod(logit_vals, 2);          
% 
%         L_zero = sum(w_mv.* gvals); 
%         log_L_zero =  log(max(L_zero, realmin)); % avoid underflow
% 
%     else % for d0>=3
%         L = chol(Sigma_cond, 'lower');
%         Z = sqrt(2) * (xk3 * L') + mu_cond'; 
%         w_mv = wk3 / (pi^(3/2));   
% 
%         logit_vals = 1 ./ (1 + exp( -(gamma0 - Z * gamma1) ));
%         gvals      = prod(logit_vals, 2);          
% 
%         L_zero = sum(w_mv.* gvals); 
%         log_L_zero = log(max(L_zero, realmin));  % avoid underflow
%     end
% 
%     logL_v = log_L_pos + log_L_zero; 
% end
