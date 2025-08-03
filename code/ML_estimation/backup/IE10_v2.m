%given that the voyage likelihood is zero I will try to get non-zero
%likelihoods by changing the parameters of the simulated data. 
clear
clc
rng(123)
% parameters 
C = 1000;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 20;     % maximum of captain voyages. 

s_omega = 2 * eye(J);     % std. dev. of product‐specific noise omega_vj

alpha  = [1.55; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.0009, 0; 0.011, 0; 0.009, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 2;                % intercept for zero‐vs‐positive
gamma1 = 2;                % slope on log(w_{1vj})


Tsim = IE10_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);

theta_real = [beta(:,1); alpha; delta; gamma0; gamma1; diag(s_omega)];

% prepare data vectors
n     = size(Tsim,1);
d     = Tsim.isPositive;
Y     = Tsim.Y_vj;
j_i   = Tsim.productID;
Xmat  = [Tsim.X1];
Tau = Tsim.Duration;
c_id = Tsim.captainID; % use built‑in captain ID

% parameter guess. 
in_beta   = [1; 1; 1];          % initial β's
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %6:8
in_gamma0 = .4;               % initial γ0 9 
in_gamma1 = 1;               % initial γ1 10
in_somega = [1; 3; .5];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];



% G = groupsummary(Tsim, 'productID', 'mean', 'isPositive')
% 
% Tpos = Tsim(Tsim.Y_vj>0, :);
% 
% G = groupsummary(Tpos, 'productID', 'mean', 'Y_vj')
% 
% r = corr(Tsim.Duration, Tsim.Y_vj);
% fprintf('Correlation(Duration,Y_vj) = %.4f\n', r);


alpha2  = [0.5; .5; .5];    % exponent 
delta2  = [3; 6; 1];    % coefficient of captain random effects
beta2   = [0.02, 0.1; 0.02, 0.1; 0.02, 0.1];        % coefficient of ship chars (weight and type).
gamma02 = 1;                % intercept for zero‐vs‐positive
gamma12 = .2;                % slope on log(w_{1vj})
Tsim2 = IE10_gen_data(C, J, Vmax, s_omega, alpha2, delta2, beta2, gamma02, gamma12);
% prepare data vectors
n2     = size(Tsim2,1);
d2     = Tsim2.isPositive;
Y2     = Tsim2.Y_vj;
j_i2   = Tsim2.productID;
Xmat2  = [Tsim2.X1, Tsim2.X2];
Tau2 = Tsim2.Duration;
c_id2= Tsim2.captainID; % use built‑in captain ID


% 
% x = repelem(10 + 0.1*randn(ceil(C*J*Vmax/3),1), 3);
% x = x(1:C*J*Vmax);
% Xmat2(:, 1) = x; 


%%

[xk, wk] = IE9_hermiteGaussRule(30);  % 30-node rule

%see backup_10 for the most voyageLik functions, we do not use them since
%captain likelihood integrates directly. 

a_c = 1.14; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 

res = voyageLik(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 



theta02 = [1; 0; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];
x_v1 = x_v * [1; 0]; 
res1 = voyageLik_prev(theta0, a_c, d_v, Y_v, x_v1, tau_v,  xk, wk); 



d_v2 = d2(1:3); 
Y_v2 = Y2(1:3); 
x_v2 = Xmat2(3,:)';
tau_v2 = Tau2(3); 
res2 = voyageLik_prev(theta02, a_c, d_v2, Y_v2, x_v2, tau_v2, xk, wk);

%%
mask  = (c_id==2);
d1    = d(mask);
Y1    = Y(mask);
X1    = Xmat(mask,:);
Tau1  = Tau(mask);

c1_lik = captainLik(theta0, d1, Y1, X1, Tau1, xk, wk);

 

%% likelihood

L_global = globalLik(theta0, d, Y, Xmat, Tau, c_id, xk, wk)


%%

% Set up negative log-likelihood for minimization
negLL = @(th) -globalLik(th, d, Y, Xmat, Tau, c_id, xk, wk);

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


