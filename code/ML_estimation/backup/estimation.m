clear
clc
rng(123)

T = readtable('../../Data/temps/clean_ML_estimation.xlsx');
disp(head(T));


% parameters 
C = numel( unique(T.captainID) );     % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = max(T.n_voy);  % maximum of captain voyages. 

% prepare data vectors
n     = size(T,1);
d     = T.d;
Y     = T.Y;
j_i   = T.productID;
Xmat  = T.X1;
Tau = T.Tau;
c_id = T.captainID; % use built‑in captain ID

% parameter guess. 
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [10 ; 10; 10]; %6:8
in_gamma0 = 1.4;               % initial γ0 9 
in_gamma1 = 20;               % initial γ1 10
in_somega = [4; 4; 4];       % initial σ_{ω,j} (σ=1)

%Pr(Y>0) approx .5 hence \gamma0 - \gamma1 z approx 0, since avg. weight is
%around 200
in_beta   = [1/2000];          % initial β's

theta0 = [in_alpha; in_delta; in_gamma0; in_gamma1; in_somega; in_beta];

%%

[xk, wk] = IE9_hermiteGaussRule(30);  % 30-node rule


function L_v = voyageLik(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk )
    J = 3;

    % Unpack parameters
    alpha       = theta(1:3);
    delta       = theta(4:6);
    gamma0      = theta(7);
    gamma1      = theta(8);
    sigma_omega = theta(9:11);
    beta = theta(12); 

    % Precompute: mu_z and sd_z for all j
    mu_z = x_v' * beta + delta * a_c;  % 3x1
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


a_c = 1.14; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:);
tau_v = Tau(3); 

res = voyageLik(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk);

%%
 
function L_c = captainLik(theta, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk)
    % Vectorized captain-level likelihood using Gauss-Hermite
    J       = 3;
    Ncap    = numel(d_cap) / J;
    ac_nodes = sqrt(2) * xk;  % nodes for a_c ~ N(0,1)

    % Reshape stacked vectors into matrices: J rows × Ncap columns
    d_mat    = reshape(d_cap, J, Ncap);
    Y_mat    = reshape(Y_cap, J, Ncap);
    tau_vec  = tau_v_cap(1:J:end)';           % 1×Ncap row vector
    K        = size(Xmat_cap,2);
    X3d      = reshape(Xmat_cap, J, Ncap, K); % J×Ncap×K

    % Unpack parameters
    alpha       = theta(1:3);
    delta       = theta(4:6);
    gamma0      = theta(7);
    gamma1      = theta(8);
    sigma_omega = theta(9:11);
    beta = theta(12); 



    % Preallocate likelihoods per quadrature node
    Mq      = length(xk);
    L_vals  = zeros(Mq, 1);

    % Loop over Gauss-Hermite nodes
    for i = 1:Mq
        a_c = ac_nodes(i);

        % Linear predictor for all voyages at once
        mu_part = zeros(J, Ncap);
        for kk = 1:K
            mu_part = mu_part + X3d(:,:,kk) * beta(kk);
        end
        mu_z_mat = mu_part + delta * a_c;            % J×Ncap
        sd_mat   = sigma_omega * ones(1, Ncap);       % J×Ncap

        % Positive-output likelihood
        u_raw = bsxfun(@rdivide, log(Y_mat), alpha);             % J×Ncap
        u_val = bsxfun(@minus, u_raw, log(tau_vec));             % J×Ncap (tau_vec is 1×Ncap)
        p_val = 1 ./ (1 + exp(gamma0 - gamma1 * u_val));
        f_z   = normpdf(u_val, mu_z_mat, sd_mat);
        jac   = 1 ./ (alpha * ones(1, Ncap) .* Y_mat);
        L_j_mat = p_val .* f_z .* jac;

        % Zero-output adjustment via GH
        zero_idx = (d_mat == 0);
        if any(zero_idx(:))
            nodes     = bsxfun(@plus, reshape(mu_z_mat, 1, J, Ncap), ...
                               reshape(sqrt(2) * sd_mat, 1, J, Ncap) .* reshape(xk, Mq, 1, 1));
            integrand = 1 ./ (1 + exp(gamma0 - gamma1 * nodes));
            weighted  = bsxfun(@times, reshape(wk, Mq, 1, 1), integrand);
            int_vals  = squeeze(sum(weighted, 1)) / sqrt(pi);        % J×Ncap
            L_j_mat(zero_idx) = 1 - int_vals(zero_idx);
        end

        % Product over products and voyages
        L_vals(i) = prod(L_j_mat(:));
    end

    % Aggregate over nodes
    L_c = sum(wk .* L_vals) / sqrt(pi);
    L_c = L_c * 1e11; 
end


% Example usage (outside this file):
% Extract for captain 1:
mask  = (c_id==11);
d1    = d(mask);
Y1    = Y(mask);
X1    = Xmat(mask,:);
Tau1  = Tau(mask);

c1_lik = captainLik(theta0, d1, Y1, X1, Tau1, xk, wk);


 

 %% likelihood
 

function L = globalLik(theta, d_vec, Y_vec, X_mat, Tau_vec, c_id, xk, wk)
% globalLik  Computes full-sample likelihood by multiplying captainLik
 
    C = max(c_id);
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
 
L_global = globalLik(theta0, d, Y, Xmat, Tau, c_id, xk, wk); 



 %%
theta0 = theta_hat; 
% Set up negative log-likelihood for minimization
negLL = @(th) -globalLik(th, d, Y, Xmat, Tau, c_id, xk, wk);

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',1000);

% Estimate theta by minimizing negative log-likelihood
[theta_hat, fval, exitflag, output] = fminunc(negLL, theta0, options);

% Final log-likelihood at estimated parameters
ll_hat = globalLik(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk);


%%


% with this initial guess 
% % parameter guess. 
% in_alpha     = [1; 1; 1];       % initial   α's (α=1)
% in_delta = [1 ; 1; 1]; %6:8
% in_gamma0 = 1.4;               % initial γ0 9 
% in_gamma1 = 1;               % initial γ1 10
% in_somega = [2; 2; 2];       % initial σ_{ω,j} (σ=1)
% in_beta   = [1/200];          % initial β's
% 
% 
% % with 100 iterations 1.10797041781158
% 1.01163479376759
% 1.00795805957937
% 1.00722636566159
% 0.992744684942064
% 1.00118268274198
% 1.41916813906242
% 0.954321522344002
% 2.02316716151763
% 2.01000400028768
% 2.00085873685703
% 0.0167134325886773
% 
% with 400 iterations 
% 1.42891346067695
% 1.00554787930375
% 1.11413967650350
% 1.01436220637498
% 0.936257488541038
% 1.01235906183003
% 1.58785807637281
% 0.649214515835945
% 2.06950767978607
% 2.05606856192170
% 1.97485840164116
% 0.0131017555823489





%%% with initial guess 
% parameter guess. 
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [10 ; 10; 10]; %6:8
in_gamma0 = 1.4;               % initial γ0 9 
in_gamma1 = 20;               % initial γ1 10
in_somega = [4; 4; 4];       % initial σ_{ω,j} (σ=1)

%Pr(Y>0) approx .5 hence \gamma0 - \gamma1 z approx 0, since avg. weight is
%around 200
in_beta   = [1/2000];   

with 300 iterations 

1.2229
0.9135
0.8881
9.9985
9.9876
9.9967
1.3997
19.9981
4.0511
4.0006
3.9883
0.0014