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

[xk, wk] = HG(15);  % 30-node rule
[xk2, wk2] = HG2D(15); 
[xk3, wk3] = HG3D(15);
%%
a_c = .76; d_v = d(1:3); Y_v = Y(1:3); x_v = Xmat(3,:)'; tau_v = Tau(3); 

a_c = .76; d_v = d(4:6); Y_v = Y(4:6); x_v = Xmat(4,:)'; tau_v = Tau(4); 

a_c = .76; d_v = d(19:21); Y_v = Y(19:21); x_v = Xmat(19,:)'; tau_v = Tau(19); 

a_c = .76; d_v = d(28:30); Y_v = Y(28:30); x_v = Xmat(28,:)'; tau_v = Tau(28); 


res = L_v_corr(theta0, a_c, d_v, Y_v, x_v, tau_v); 
res2 = L_v_corr_v2(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3); 

theta2 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; diag(in_somega)];
res3 = voyageLik(theta2, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 

%%

a_c = 76; 
d_c = d(1: J*Vmax); 
Y_c = Y(1: J*Vmax);
x_c = Xmat(1: J*Vmax);
tau_c = Tau(1: J*Vmax); 
theta = theta0; 

%%
 

res5 = L_c_corr(theta, a_c, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3)

tic
res6 = L_c_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc


tic
rest7 = captainLik_full(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
toc


res8 = globalLik_corr(theta, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 