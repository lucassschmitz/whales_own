%given that the voyage likelihood is zero I will try to get non-zero
%likelihoods by changing the parameters of the simulated data. 
clear
clc
rng(123)
% parameters 
C = 3000;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 20;     % maximum of captain voyages. 

s_omega = 2 * eye(J);     % std. dev. of product‐specific noise omega_vj
s_omega = [1,1,1; 1,1,1; 1,1,1];

alpha  = [1.55; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.0009, 0; 0.011, 0; 0.009, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 2;                % intercept for zero‐vs‐positive
gamma1 = 2;                % slope on log(w_{1vj})


Tsim = IE10_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);

theta_real = [beta(:,1); alpha; delta; gamma0; gamma1; s_omega(:)];

% check the moments 
% G = groupsummary(Tsim, 'productID', 'mean', 'isPositive')
% Tpos = Tsim(Tsim.Y_vj>0, :);
% G = groupsummary(Tpos, 'productID', 'mean', 'Y_vj')
% r = corr(Tsim.Duration, Tsim.Y_vj);
% fprintf('Correlation(Duration,Y_vj) = %.4f\n', r);
% 

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
in_somega = [1, 0, 0; 0, 3, 0; 0, 0, .5];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:)];

[xk, wk] = HG(30);  % 30-node rule
[nodes, weights] = HG3D(4);
%%
clc
a_c = .76; 
d_v = d(4:6); 
Y_v = Y(4:6); 
x_v = Xmat(4,:)';
tau_v = Tau(4); 

W = 2 + randn(3,1); 

res = voyageLik_initial0(theta0, d_v, W, x_v);
res2 = voyageLik_initial(theta0, d_v, W, x_v); 

aux = l_v_integration(theta0, d_v, W, x_v); 
%%

function L_v = voyageLik_initial0(theta, d_v, W, x_v) 
    % allows 1. product specific betas and 2. correlation across products     

    J = 3;
    p = numel(x_v);

    idx      = J*p;

     % Unpack remaining parameters
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);
    
    L_v = 1; % initialize the likelihood
    for j = 1: J 
        aux1 = exp(gamma0 - gamma1*W(j))^(1-d_v(j)); 
        aux2 = 1 + exp(gamma0 - gamma1 * W(j)); 
        L_v  = L_v * (aux1/aux2);
    end
end 

function L_v = voyageLik_initial(theta, d_v, w_hat_v, x_v)

J = numel(d_v);
p = numel(x_v);

%—— unpack γ0 and γ1 from theta
idx     = J*p + 2*J;
gamma0  = theta(idx + 1);
gamma1  = theta(idx + 2);

%—— latent index η_j = γ0 – γ1·ŵ_{vj}
eta = gamma0 - gamma1 .* w_hat_v;    % J×1

%—— per‐product likelihood factors
num = exp(eta) .^ (1 - d_v);  % numerator = exp(η)^1−d_j
den = 1 + exp(eta);           % denominator = 1+exp(η)
L_v = prod(num ./ den);
end

function [L_v, f_joint] = l_v_integration(theta, a_c, d_v, w_hat_v, x_v) 
    % allows 1. product specific betas and 2. correlation across products     
    J = 3;
    p = numel(x_v);

    beta_vec = theta(1:J*p);
    beta_mat = reshape(beta_vec, J, p);      % J×p matrix of betas
    idx      = J*p;

    % Unpack remaining parameters
    alpha       = theta(idx+1    : idx+J);         % J×1          
    delta       = theta(idx+J+1  : idx+2*J);       % J×1          
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);

    Omega_vec = theta(idx+2*J+3 : idx+2*J+2+J^2);
    sigma_omega   = reshape(Omega_vec, J, J);
    
    % Precompute: mu_z and sd_z for all j
    mu_z = beta_mat *x_v + delta * a_c;  % 3x1
    sd_z = sigma_omega;

    f_joint = mvnpdf(u_val', mu_z', sigma_omega) 





end 







function L_v = voyageLik_v2(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk )
    % allows 1. product specific betas and 2. correlation across products     
    J = 3;
    p = numel(x_v);

    beta_vec = theta(1:J*p);
    beta_mat = reshape(beta_vec, J, p);      % J×p matrix of betas
    idx      = J*p;

     % Unpack remaining parameters
    alpha       = theta(idx+1    : idx+J);         % J×1          
    delta       = theta(idx+J+1  : idx+2*J);       % J×1          
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);

    Omega_vec = theta(idx+2*J+3 : idx+2*J+2+J^2);
    sigma_omega   = reshape(Omega_vec, J, J);
    
   



    % Precompute: mu_z and sd_z for all j
    mu_z = beta_mat *x_v + delta * a_c;  % 3x1
    sd_z = sigma_omega;

    % Precompute common terms
    u_val  = log(Y_v) ./ alpha - log(tau_v);
    p_val  = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));

    f_joint = mvnpdf(u_val', mu_z', sigma_omega); 
    %f_z    = normpdf(u_val, mu_z, sd_z);
    jac    = 1 ./ (alpha .* Y_v);

    % Vectorized likelihoods for d_v == 1 and d_v == 0
    L_j = zeros(J,1);

    pos_idx = d_v == 1;
    L_j(pos_idx) = p_val(pos_idx) .* f_z(pos_idx) .* jac(pos_idx);

    % Use Gauss-Hermite quadrature for the zero-output case
    for j = find(~pos_idx)'
        % Change of variables: u = sqrt(2)*sd*xk + mu
        nodes = sqrt(2) * sd_z(j) * xk + mu_z(j);
        weights = wk .* (1 ./ (1 + exp(gamma0 - gamma1 * nodes)));
        L_j(j) = 1 - sum(weights) / sqrt(pi);
    end

    % Final voyage likelihood is product across j
    L_v = prod(L_j);
end



%%
a_c = .76; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 
 
res = voyageLik(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 

%%
mask  = (c_id==284);
d_cap    = d(mask);
Y_cap    = Y(mask);
Xmat_cap    = Xmat(mask,:);
tau_v_cap  = Tau(mask);

LogLc = captainLik_v2(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk);        

%%
L_global3 = globalLik_v3(theta0, d, Y, Xmat, Tau, c_id, xk, wk);
%%

% Set up negative log-likelihood for minimization
negLL = @(th) -globalLik_v3(th, d, Y, Xmat, Tau, c_id, xk, wk);

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',4000, ...
    'MaxFunctionEvaluations',4000);

% Estimate theta by minimizing negative log-likelihood
[theta_hat, fval, exitflag, output] = fminunc(negLL, theta0, options);

% Final log-likelihood at estimated parameters
ll_hat = globalLik(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk);

dist = norm(theta_hat - theta_real);

%%

% Lower bounds for non-negativity on alpha (3:5),  delta (6:8), gamma1(10), sigma_omega(11:13)
lb = -inf(size(theta0));
lb([3:5, 6:8,10:13]) = 0;

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'TolFun',1e-4, 'TolX',1e-4, 'OptimalityTolerance',1e-4, ...
    'MaxIterations',300, 'MaxFunctionEvaluations',300);

[theta_hat2, fval, exitflag, output] = fmincon(negLL, theta0, [], [], [], [], lb, [], [], options);
dist2 = norm(theta_hat2 - theta_real);





function [x, w] = HG(n)
    i = (1:n-1)';
    a = sqrt(i/2);
    T = diag(a,1) + diag(a,-1);
    [V, D] = eig(T);
    [x_sorted, idx] = sort(diag(D));     % GH nodes
    V = V(:, idx);
    x = x_sorted;
    w = (sqrt(pi) * (V(1,:).^2))';       % ensure w is column (n×1)
end


function [nodes, weights] = HG3D(n)
    %returns an (n^3×3) array of 3D nodes and an (n^3×1) vector of weights

    [x, w] = HG(n);

    % Build 3‐D grids of nodes and weights
    [X1, X2, X3] = ndgrid(x, x, x);
    [W1, W2, W3] = ndgrid(w, w, w);

    % Pack into list of 3D nodes and corresponding weights
    nodes   = [X1(:), X2(:), X3(:)];
    weights = W1(:) .* W2(:) .* W3(:);
end



function logLc = captainLik_v2(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk)
    % CAPTAIN-LEVEL likelihood using Gauss-Hermite with product-specific betas
    % and using logs to avoid underflowing. Returns the log captain
    % likelihood .
    J       = 3;
    Ncap    = numel(d_cap) / J;
    ac_nodes = sqrt(2) * xk;  % nodes for a_c ~ N(0,1)

    % Reshape stacked vectors
    d_mat    = reshape(d_cap, J, Ncap);
    Y_mat    = reshape(Y_cap, J, Ncap);
    tau_vec  = tau_v_cap(1:J:end)';             % 1×Ncap
    p        = size(Xmat_cap,2);
    X3d      = reshape(Xmat_cap, J, Ncap, p);   % J×Ncap×p

    % Unpack parameters
    beta_vec    = theta(1 : J*p);              
    beta_mat    = reshape(beta_vec, J, p);      % J×p matrix  
    idx         = J*p;                          

    alpha       = theta(idx+1    : idx+J);
    delta       = theta(idx+J+1  : idx+2*J);
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);
    sigma_omega = theta(idx+2*J+3: idx+3*J+2);

    Mq      = length(xk);
    logL_vals = zeros(Mq,1);             %stores logs to avoid underflow. 

    % Loop over Gauss-Hermite nodes
    for i = 1:Mq
        a_c = ac_nodes(i);

        % compute means for each product & voyage
        mu_part = zeros(J, Ncap);
        for kk = 1:p
            mu_part = mu_part + bsxfun(@times, X3d(:,:,kk), beta_mat(:,kk));
        end
        mu_z_mat = mu_part + delta * a_c;          % J×Ncap
        sd_mat   = sigma_omega * ones(1, Ncap);

        % per-product×voyage likelihoods
        u_raw   = bsxfun(@rdivide, log(Y_mat), alpha);
        u_val   = bsxfun(@minus, u_raw, log(tau_vec));
        p_val   = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
        f_z     = normpdf(u_val, mu_z_mat, sd_mat);
        jac     = 1 ./ (alpha * ones(1,Ncap) .* Y_mat);
        L_j_mat = p_val .* f_z .* jac;

        % zero-output cells replacement
        zero_idx = (d_mat == 0);
        if any(zero_idx(:))
            nodes     = bsxfun(@plus, reshape(mu_z_mat,1,J,Ncap), ...
                               reshape(sqrt(2)*sd_mat,1,J,Ncap).*reshape(xk,Mq,1,1));
            integrand = 1 ./ (1 + exp(gamma0 - gamma1 * nodes));
            weighted  = bsxfun(@times, reshape(wk,Mq,1,1), integrand);
            int_vals  = squeeze(sum(weighted,1)) / sqrt(pi);
            L_j_mat(zero_idx) = 1 - int_vals(zero_idx);
        end

        %L_j_mat = max( L_j_mat, realmin ); % floor negative values. 

        % accumulate log-likelihood
        logL_vals(i) = sum(log(L_j_mat(:)));     
    end
   
    % integrate via log-sum-exp
    log_w = log(wk(:)) - 0.5*log(pi);            
    log_tk = log_w + logL_vals;                 
    M      = max(log_tk);                      
    sum_exp = sum(exp(log_tk - M));            
    logLc  = M + log(sum_exp);                
    %L_c    = exp(logLc);                       
end

function Tsim = IE10_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1)

%differences with IE9_gen_data: 
% 1.  we tweak the function to produce data similar to the real data
% 2. we allow for product specific beta. 
% 3. for simplicity duration was kept fixed to 1, now I allow for it to
% change production 

%   beta      - J×2 matrix of covariate coefficients (product-specific)  <-- CHANGED


% parameters obtained from our data 
%1. duration (taken as continuous; mean 2.5, s.d 1.5, min 1) 
mu_tau = 2.5; 
sd_tau = 1.5; 
min_tau = 1; 

%2. tonnage (mean 272, sd 103, min 45 )
mu_x = 272;
sd_x   = 103;
min_x  = 15; % set 15 to get around 45 



% Pre‐allocate the “long” vectors exactly as before
N  = C * Vmax * J;
captainID = zeros(N,1);
Nvoy = zeros(N,1);
productID = zeros(N,1);

X1        = zeros(N,1);
X2        = zeros(N,1);

a_c       = zeros(N,1);
omega_vj  = zeros(N,1);   % STILL a scalar per row, but filled from a J‐vector
log_wvj   = zeros(N,1);
P_pos     = zeros(N,1);
isPositive= zeros(N,1);
Y_vj      = zeros(N,1);
Dur       = zeros(N,1); 

idx = 0;


% Pre‐draw all a_c for each captain
true_ac = randn(C,1);

for c = 1:C
    for v = 1:Vmax
        

        % Draw a J‐vector of omegas for this (captain c, voyage v):
        Omega_vec = mvnrnd( zeros(1,J), s_omega );  
        
        % voyage duration 
        mu_Z = mu_tau - min_tau; 
        sigma2 = log(1 + (sd_tau/mu_Z)^2); 
        mu = log(mu_Z) - sigma2/2; 
        Dur_aux = lognrnd(mu, sqrt(sigma2), 1, 1); 
        Dur_aux = Dur_aux+ min_tau;

        idx_start = J * (Vmax * (c-1) + (v-1)) + 1;
        idx_end = idx_start + J - 1;

        Dur(idx_start:idx_end) = Dur_aux;

        % ship chars. 
        mu_Z = mu_x - min_x; 
        sigma2 = log(1 + (sd_x/mu_Z)^2); 
        mu = log(mu_Z) - sigma2/2; 
        x1_v = lognrnd(mu, sqrt(sigma2), 1, 1); 
        x1_v = x1_v + min_x; 

        x2_v = rand > 0.5; 

        % 2) Loop over products j = 1..J, assign each component of Omega_vec:
        for j = 1:J
            idx = idx + 1;

            captainID(idx) = c;
            Nvoy(idx)  = v;
            productID(idx)  = j;

            % % Covariates X1, X2
            X1(idx) = x1_v;
            X2(idx) = x2_v;

            % Captain skill, replicated for every row with the same captain c:
            a_c(idx) = true_ac(c);

            % Now assign the j‐th component of the MVN draw:
            omega_vj(idx) = Omega_vec(j);

            % Compute log w_{vj} = X_f * beta + delta_j * a_c + omega_vj:
            log_wvj(idx) = beta(j,1)*x1_v + beta(j,2)*x2_v + delta(j)*a_c(idx) + omega_vj(idx);
                          
            P_pos(idx) = 1 / (1 + exp(gamma0 - gamma1 * log_wvj(idx)));
            
            isPositive(idx) = (rand < P_pos(idx));

            tau_v = Dur(idx);  
            if isPositive(idx)
                Y_vj(idx) = (tau_v * exp(log_wvj(idx)))^alpha(j);
            else
                Y_vj(idx) = 0;
            end
        end
    end
end

 
% Finally, pack into a table as before:
Tsim = table( ...
   captainID, Nvoy, productID, ...
   X1, X2, ...
   a_c, omega_vj, log_wvj, P_pos, isPositive, Y_vj, Dur, ...
   'VariableNames', { ...
     'captainID','Nvoy','productID', ...
     'X1','X2', ...
     'a_c','omega_vj','log_wvj','P_pos','isPositive','Y_vj', 'Duration' } );


end

function L = globalLik_v3(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk)
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
        logL_c(c) = captainLik_v2(theta, d_c, Y_c, X_c, Tau_c, xk, wk);
    end
    
    % sum up all the per‐captain logs
    L = sum(logL_c);
end

function L_v = voyageLik_prev(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk )
    % allows 1. product specific betas
    J = 3;
    p = numel(x_v);

    beta_vec = theta(1:J*p);
    beta_mat = reshape(beta_vec, J, p);      % J×p matrix of betas
    idx      = J*p;

     % Unpack remaining parameters
    alpha       = theta(idx+1    : idx+J);         % J×1          
    delta       = theta(idx+J+1  : idx+2*J);       % J×1          
    gamma0      = theta(idx+2*J+1);
    gamma1      = theta(idx+2*J+2);
    sigma_omega = theta(idx+2*J+3: idx+3*J+2);      % J×1         



    % Precompute: mu_z and sd_z for all j
    mu_z = beta_mat *x_v + delta * a_c;  % 3x1
    sd_z = sigma_omega;

    % Precompute common terms
    u_val  = log(Y_v) ./ alpha - log(tau_v);
    p_val  = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
    f_z    = normpdf(u_val, mu_z, sd_z);
    jac    = 1 ./ (alpha .* Y_v);

    % Vectorized likelihoods for d_v == 1 and d_v == 0
    L_j = zeros(J,1);

    pos_idx = d_v == 1;
    L_j(pos_idx) = p_val(pos_idx) .* f_z(pos_idx) .* jac(pos_idx);

    % Use Gauss-Hermite quadrature for the zero-output case
    for j = find(~pos_idx)'
        % Change of variables: u = sqrt(2)*sd*xk + mu
        nodes = sqrt(2) * sd_z(j) * xk + mu_z(j);
        weights = wk .* (1 ./ (1 + exp(gamma0 - gamma1 * nodes)));
        L_j(j) = 1 - sum(weights) / sqrt(pi);
    end

    % Final voyage likelihood is product across j
    L_v = prod(L_j);
end