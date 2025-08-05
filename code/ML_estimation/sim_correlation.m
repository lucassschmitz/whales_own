clear
clc

rng(123)
% parameters 
C = 50;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 25;     % maximum of captain voyages. 

s_omega = 2 * eye(J);     % std. dev. of product‐specific noise omega_vj
s_omega = [2 1 1; 1 2 1; 1 1 2];
disp(['Lowest eigenvalue = ' num2str(min(eig(s_omega)))]);

 
alpha  = [1.3; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.0019, 0; 0.021, 0; 0.029, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 4;                % intercept for zero‐vs‐positive
gamma1 = 2;                % slope on log(w_{1vj})


Tsim = IE10_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);

theta_real = [beta(:,1); alpha; delta; gamma0; gamma1; s_omega(:)];


% prepare data vectors
n     = size(Tsim,1);
d     = Tsim.isPositive;
Y     = Tsim.Y_vj;
j_i   = Tsim.productID;
Xmat  = [Tsim.X1];
Tau = Tsim.Duration;
c_id = Tsim.captainID; % use built‑in captain ID

% parameter guess. 
in_beta   = [0.01; 0.01; 0.01];          % initial β's
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %6:8
in_gamma0 = .4;               % initial γ0 9 
in_gamma1 = 1;               % initial γ1 10
in_somega = [1, 0.5, 0; 0.5, 3, 0; 0, 0, .5];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:)];
theta_chol0 = to_chol_theta(theta0); % theta reduced 

[xk, wk] = HG(15);  % 30-node rule
[xk2, wk2] = HG2D(15); 
[xk3, wk3] = HG3D(15);

% CODE USED TO CHECK THE FUNCTIONS 
a_c = .76; d_v = d(1:3); Y_v = Y(1:3); x_v = Xmat(3,:)'; tau_v = Tau(3); 
a_c = .76; d_v = d(4:6); Y_v = Y(4:6); x_v = Xmat(4,:)'; tau_v = Tau(4); 
a_c = .76; d_v = d(19:21); Y_v = Y(19:21); x_v = Xmat(19,:)'; tau_v = Tau(19); 
a_c = .76; d_v = d(28:30); Y_v = Y(28:30); x_v = Xmat(28,:)'; tau_v = Tau(28); 
a_c = .76; d_v = d(46:48); Y_v = Y(46:48); x_v = Xmat(46,:)'; tau_v = Tau(28); 
a_c = .76; d_v = d(157:159); Y_v = Y(157:159); x_v = Xmat(157,:)'; tau_v = Tau(157); 

% tic
% res = L_v_corr(theta0, a_c, d_v, Y_v, x_v, tau_v); 
% toc
tic
res2 = L_v_corr_v2(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3); 
toc
tic
res21 = exp(L_v_corr_v3(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3) ); 
toc
%% 

theta2 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; diag(in_somega)];
res3 = voyageLik(theta2, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 

mask   = (c_id == 3);
d_c    = d(mask);
Y_c    = Y(mask);
x_c    = Xmat(mask, :);
tau_c  = Tau(mask);
a_c = 0; 
theta = theta0; 
tic
res5 = L_c_corr(theta0, a_c, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc
%%
tic
res6 = L_c_corr_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc
tic
res7 = L_c_corr_int_v2(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc
tic
res71 = L_c_corr_int_v3(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc
abs(res6- res7) < 1e-10    & abs(res6- res71) < 1e-10 

tic
res8 = globalLik_corr(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
toc
tic 
[res9, B] = globalLik_corr_v2(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
toc
abs(res8- res9) < 1e-10


%% Estimation 
clc
negLL_red = @(chol_theta) -globalLik_corr_v2(to_theta(chol_theta), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

% % Diagnostic1
% % Add this after creating theta_chol0
% fprintf('--- Starting Diagnostic ---\n');
% initial_negLL = negLL_red(theta_chol0);
% 
% fprintf('Initial objective function value: %f\n', initial_negLL);
% if isnan(initial_negLL) || isinf(initial_negLL)
%     fprintf('WARNING: Objective function returned NaN or Inf!\n');
% end
% fprintf('--- End Diagnostic ---\n');
% 
% theta_chol02 = theta_chol0; 
% theta_chol02(1:3) = [.2; .2; .2];
% initial_negLL2 = negLL_red(theta_chol02);
%
% %chekc that objective is not flat 
% a = negLL_red(theta_chol0) 
% b = negLL_red(theta_chol0+1e-4)


% % --- Diagnostic2: Plot Likelihood Surface ---
% param_to_test = 1; % Test the first parameter
% theta_start = theta_chol0;
% n_points = 50;
% param_range = linspace(theta_start(param_to_test) * 0.9, theta_start(param_to_test) * 1.1, n_points);
% y_vals = zeros(n_points, 1);
% 
% for i = 1:n_points
%     theta_test = theta_start;
%     theta_test(param_to_test) = param_range(i);
%     y_vals(i) = negLL_red(theta_test);
% end
% 
% figure;
% plot(param_range, y_vals);
% title(['Log-Likelihood vs. Parameter ', num2str(param_to_test)]);
% xlabel(['Value of Parameter ', num2str(param_to_test)]);
% ylabel('Negative Log-Likelihood');

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',200, ... 
    'MaxFunctionEvaluations', 30, ...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',      1e-15) 
 
[theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );

theta_hat = to_theta(theta_red_hat); % recover the full vector 
ll_hat    = globalLik_corr(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); % evaluate function. 

%save('est_corr_unrestricted', 'theta_hat')

%% Initial values
% Number of different starting points
num_inits = 5;

% Store results in a table for easier comparison
results_table = table('Size', [num_inits, 4], 'VariableTypes', {'double', 'cell', 'cell', 'double'}, ...
                      'VariableNames', {'InitialValueID', 'EstimatedTheta', 'DifferenceFromReal', 'Distance'});

% Generate different initial values by adding noise directly to the Cholesky parameters
theta_chol_inits = repmat(theta_chol0, 1, num_inits) + randn(size(theta_chol0, 1), num_inits) * 0.2;

%% 
for i = 1:num_inits
    fprintf('Running Minimization for Initial Value %d...\n', i); 


    aux = globalLik_corr(to_theta(theta_chol_inits(:, i)), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3)

    % Run the optimizer
    [theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol_inits(:, i), options );
    
    % Store results
    theta_hat = to_theta(theta_red_hat);
    difference = theta_hat - theta_real;
    distance = norm(difference);
    
    results_table.InitialValueID(i) = i;
    results_table.EstimatedTheta{i} = theta_hat';
    results_table.DifferenceFromReal{i} = difference';
    results_table.Distance(i) = distance;
end

%% Comparison of Estimates with theta_real

fprintf('\n\n--- Comparison of Final Estimates with True Theta ---\n\n');
fprintf('Theta Real:\n');
disp(theta_real');
disp('--- Estimation Results ---');
disp(results_table);


% 
% 
% function L = globalLik_corr_v2(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk, xk2, wk2, xk3, wk3)
% % globalLik_corr  Full-sample log-likelihood across all captains allowing
% % for correlation. 
% %   theta    parameter vector as above
% %   d_vec,Y_vec,X_mat,Tau_vec  data stacks [J·Ntotal×1 or ×p]
% %   c_id     [J·Ntotal×1] captain IDs (1...C)
% %   xk,wk,...  GH nodes & weights
%     %disp(theta(end-9:end)')
%     C = max(c_id);
%     logL_c = zeros(C,1);
% 
%     % ——— START parallel pool if possible ———
%     try
%         poolobj = gcp;
%         if isempty(poolobj)
%             parpool;
%         end
%     catch
%         if exist('matlabpool','file') && matlabpool('size')==0
%             matlabpool open;
%         end
%     end
% 
%     % Compute each captain's log-likelihood in parallel
%     parfor c = 1:C
%         mask   = (c_id == c);
%         d_c    = d_vec(mask);
%         Y_c    = Y_vec(mask);
%         X_c    = X_mat(mask, :);
%         Tau_c  = Tau_vec(mask);
% 
%         logL_c(c) = L_c_corr_int_v3(theta, d_c, Y_c, X_c, Tau_c, ...
%                                     xk, wk, xk2, wk2, xk3, wk3);
%     end
% 
%     % Sum across captains
%     L = sum(logL_c);
% end
% 
% 
% 
% 
% 
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
% 
% 
% function logL_v  = L_v_corr_v3(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3)
%     % implements L_v_corr (produces same result up to the log function) but 
%     % trying to optimize the integration, spetially when the voyage produces 0
%     % output the 3D integral is hard to compute. 
%     % implements log voyage likelihood given with correlation  
% 
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
% 
% 
% 
% 
