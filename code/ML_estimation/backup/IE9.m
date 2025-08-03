clear
clc
rng(123)
% parameters 
C = 500;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 20;     % maximum of captain voyages. 

s_omega = 2 * eye(J);     % std. dev. of product‐specific noise omega_vj

alpha  = [0.8; .7; .9];    % exponent 
delta  = [3; 6; 1];    % coefficient of captain random effects
beta   = [0.3; 0.1];        % coefficient of ship chars (weight and type).
gamma0 = 1.4;                % intercept for zero‐vs‐positive
gamma1 = .2;                % slope on log(w_{1vj})


Tsim = IE9_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);
theta_real = [beta; alpha; delta; gamma0; gamma1; diag(s_omega)];

% prepare data vectors
n     = size(Tsim,1);
d     = Tsim.isPositive;
Y     = Tsim.Y_vj;
j_i   = Tsim.productID;
Xmat  = [Tsim.X1, Tsim.X2];
Tau = Tsim.Duration;
c_id = Tsim.captainID; % use built‑in captain ID

% parameter guess. 
in_beta   = [1; 1];          % initial β's
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %6:8
in_gamma0 = 1.4;               % initial γ0 9 
in_gamma1 = 1;               % initial γ1 10
in_somega = [2; 2; 2];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];



%%

[xk, wk] = IE9_hermiteGaussRule(30);  % 30-node rule


a_c = 1.14; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 

res = voyageLik(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 

%%
 
% Example usage (outside this file):
% Extract for captain 1:
mask  = (c_id==3);
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

%% 
dist = norm(theta_hat - theta_real);

% when real params are : 
% alpha  = [0.8; 1.1; 1.0];    % exponent 
% delta  = [1.5; 1.05; 0.9];    % coefficient of captain random effects
% beta   = [0.3; 0.1];        % coefficient of ship chars (weight and type).
% gamma0 = 1.4;                % intercept for zero‐vs‐positive
% gamma1 = 1;                % slope on log(w_{1vj})
% 

% with C = 10 and max_iterations/functions 400 -> dist = 7.4 
% with C = 30 and max_iterations/functions 400 -> dist =  5.33
% with C = 30 and max_iterations/functions 200 -> dist =  1.9 
% with C = 100 and max_iterations/functions 200 -> dist =  1.9804
% with C = 300 and max_iterations/functions 200 -> dist =  1.794
% with C = 500, Vmax = 20  and max_iterations/functions 300 -> dist = 3.1071


%%

% Lower bounds for non-negativity on alpha (3:5),  delta (6:8), gamma1(10), sigma_omega(11:13)
lb = -inf(size(theta0));
lb([3:5, 6:8,10:13]) = 0;


% Optimization options using fmincon for bounds
options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'TolFun',1e-4, 'TolX',1e-4, 'OptimalityTolerance',1e-4, ...
    'MaxIterations',300, 'MaxFunctionEvaluations',300);

% Run constrained minimization
[theta_hat2, fval, exitflag, output] = fmincon(negLL, theta0, [], [], [], [], lb, [], [], options);
dist2 = norm(theta_hat2 - theta_real);

% with C = 35 and max_iterations/functions 300 -> dist = 5.1138
% with C = 500 and max_iterations/function 300 -> dist = 2.5