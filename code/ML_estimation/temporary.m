clear
clc

rng(123)
% parameters 
C = 40;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 20;     % maximum of captain voyages. 

s_omega = [2 1 1; 1 2 1; 1 1 2];
disp(['Lowest eigenvalue = ' num2str(min(eig(s_omega)))]);
 
alpha  = [1.55; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.0019, 0; 0.021, 0; 0.029, 0 ];        % coefficient of ship chars (weight and type).
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
in_beta   = [0.1; 0.1; 0.1];          % initial β's
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

%% Estimation 

negLL_red = @(chol_theta) -globalLik_corr(to_theta(chol_theta), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',100, ... 
    'MaxFunctionEvaluations', 100, ...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',      1e-10, ...
    'FiniteDifferenceStepSize',1e-4) %bigger step size to avoid zero gradietn. 
 
[theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );


theta_hat = to_theta(theta_red_hat); % recover the full vector 
ll_hat    = globalLik_corr(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); % evaluate function. 

%%
opts2 = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8);
[theta_red_hat2, fval2] = fminsearch(negLL_red, theta_chol0, opts2);

theta_red_hat2_stored = [-13.5524305275554; -6.49782773370314; -5.00933102558122; ...
    19.7786073849688;-7.15979192729329;6.29699612733020;33.8807761277274; -6.21347393762378; ... 
    88.0086625422573;10.9602451339311;0.0722344418592385;-23.5456095404422;100.933926038650; ... 
    21.9963457950554;-7.48164384085295; -0.0587296136334158;0.107493293078392]

%% 


[theta_red_hat2, fval2, exitflag2, output2] = fminunc(negLL_red, theta_red_hat2_stored, options );
