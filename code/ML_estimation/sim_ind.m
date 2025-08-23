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

alpha  = [1.55; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.0009, 0; 0.011, 0; 0.009, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 2;                % intercept for zero‐vs‐positive
gamma1 = 2;                % slope on log(w_{1vj})


Tsim = simulate_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);

theta_real = [beta(:,1); alpha; delta; gamma0; gamma1; diag(s_omega)];

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
in_somega = [1; 3; .5];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];

[xk, wk] = HG(30);  % 30-node rule

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

LogLc = captainLik(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk);        

%%
L_global3 = globalLik(theta0, d, Y, Xmat, Tau, c_id, xk, wk);
%%

% Set up negative log-likelihood for minimization
negLL = @(th) -globalLik(th, d, Y, Xmat, Tau, c_id, xk, wk);

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',4000, ...
    'MaxFunctionEvaluations',4);

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
    'MaxIterations',300, 'MaxFunctionEvaluations',2);

[theta_hat2, fval, exitflag, output] = fmincon(negLL, theta0, [], [], [], [], lb, [], [], options);
dist2 = norm(theta_hat2 - theta_real);
