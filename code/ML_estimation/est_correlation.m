 %   USE THE CODE OF SIM_CORR TO ESTIMATE THE MODEL. 
%%% Code to estimate the parameters of the model using the real data. 

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
in_beta   = [0.01; 0.01; 0.01];          % initial β's
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %7:9
in_gamma0 = .4;               % initial γ0 10 
in_gamma1 = 1;               % initial γ1 11
in_somega = [1, 0.5, .5; .5, 3, .5; .5, .5, .5];     
% σ_a is normalized to 1 

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:)];

[xk, wk] = HG(30);  % 30-node rule
[xk2, wk2] = HG2D(15); 
[xk3, wk3] = HG3D(15);

a_c = .76; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 

res = L_v_corr_v2(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);

%% 



mask  = (c_id==284);
d_cap    = d(mask);
Y_cap    = Y(mask);
Xmat_cap    = Xmat(mask,:);
tau_v_cap  = Tau(mask);

LogLc = captainLik_full(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk, xk2, wk2, xk3, wk3);        

global_lik = globalLik_corr(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 

%%

% when minimizing we need to include the restriction that the covariance
% matrix is symmetric. denote theta_red the vector without the lower
% triangle. 

theta_red0 = from_theta(theta0); 

negLL_red = @(theta_red) -globalLik_corr(to_theta(theta_red), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',20, ...
    'MaxFunctionEvaluations',20);


[theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_red0, options );

theta_hat = to_theta(theta_red_hat); % recover the full vector 
ll_hat    = globalLik_corr(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); % evaluate function. 


save('est_corr__unrestricted', 'theta_hat')


%% Loop for different initial guess. 

theta_0aux = theta0;  
theta_0aux(4) = 1.1; 
theta0_mat(:,1) = theta_0aux; 


theta_0aux = theta0 ; 
theta_0aux(5) = 1.1; 
theta_0aux(13) = 2;  
theta0_mat(:,2) = theta_0aux; 


% initialize
N = size(theta0_mat, 2);                 % number of iterations

A(N) = struct();     % preallocate struct array

for i = 1:N
    theta0 = theta0_mat(:, i);

    [theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, from_theta(theta0), options);

    % Store results
    A(i).theta_hat = theta_red_hat;   % column vector
    A(i).fval      = fval;        % scalar
    A(i).exitflag  = exitflag;    % scalar
    A(i).output    = output;      % struct (e.g., iterations, message, etc.)
end

save('est_corr_unrestrictedMat.mat', 'A');

%%
theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:)];

theta0_red = from_theta(theta0);


% Lower bounds on theta_reduced for non-negativity on alpha (4:6),  delta (7:9), gamma1(11), sigma_omega(12:14)
lb = -inf(size(theta_red0));
lb([4:6, 7:9,11, 12:14]) = 0;

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'TolFun',1e-4, 'TolX',1e-4, 'OptimalityTolerance',1e-4, ...
    'MaxIterations',10, 'MaxFunctionEvaluations',3000);

[theta_red_hat2, fval, exitflag, output] = fmincon(negLL_red, theta0_red, [], [], [], [], lb, [], [], options);


save('est_corr__restricted', 'theta_hat2')


A2(N) = struct();     % preallocate struct array

for i = 1:N
    theta0 = theta0_mat(:, i);

    [theta_red_hat2, fval, exitflag, output] = fmincon(negLL_red, theta0, [], [], [], [], lb, [], [], options);

    % Store results
    A2(i).theta_hat = theta_hat2;   % column vector
    A2(i).fval      = fval;        % scalar
    A2(i).exitflag  = exitflag;    % scalar
    A2(i).output    = output;      % struct (e.g., iterations, message, etc.)
end

save('est_corr_restrictedMat.mat', 'A2');


%% minimization stability 
 

%── 1) Load saved estimates ──────────────────────────────────────────────
load('est_corr_unrestricted.mat','theta_hat')    % fminunc initial
load('est_corr_restricted.mat','theta_hat2')  % fmincon initial
load('est_corr_unrestrictedMat.mat','A')    % struct array from fminunc runs
load('est_corr_restrictedMat.mat','A2')  % struct array from fmincon runs

%── 2) Stack into one big matrix ────────────────────────────────────────
N = numel(A);               % number of additional inits per optimizer
p = numel(theta_hat);       % dimension of theta

M = 2*(N+1);                % total number of runs
allTheta = zeros(M, p);

% fminunc block
allTheta(1,   :) = theta_hat';
for i = 1:N
    allTheta(1+i, :) = A(i).theta_hat';
end

% fmincon block
base = N+1;
allTheta(base+1,   :) = theta_hat2';
for i = 1:N
    allTheta(base+1+i, :) = A2(i).theta_hat';
end

%── 3) Compute dispersion statistics ─────────────────────────────────────
mean_theta   = mean(allTheta, 1);       % 1×p
std_theta    = std(allTheta, 0, 1);     % 1×p
var_theta    = var(allTheta, 0, 1);     % 1×p
cov_theta    = cov(allTheta);           % p×p covariance matrix
cv_theta     = std_theta ./ abs(mean_theta);      % 1×p coefficient of variation

%── 4) Display a summary table ──────────────────────────────────────────
paramNames = { ...
    'beta_bone', 'beta_sperm', 'beta_oil', ...     % β parameters
    'alpha_bone','alpha_sperm','alpha_oil', ...    % α parameters
    'delta_bone','delta_sperm','delta_oil', ...    % δ parameters
    'gamma0','gamma1', ...                         % γ parameters
    'sigma_bone','sigma_sperm','sigma_oil' ...    % σ_ω parameters
}'; %paramNames = compose("theta%02d", 1:p)';  
params_tab = table(mean_theta', std_theta', var_theta', cv_theta', ...
          'VariableNames', {'Mean','StdDev','Variance','CV'}, ...
          'RowNames', paramNames);
disp(params_tab)


%% model fit. 

% moments in the data 
% assume T is your table
productIDs = unique(T.productID);


nProd      = numel(productIDs);

% pre‐allocate
meanY   = zeros(nProd,1);
shareD1 = zeros(nProd,1);

for i = 1:nProd
    pid      = productIDs(i);
    mask     = T.productID == pid;
    meanY(i) = mean( T.Y(mask) );
    % if d is already 0/1 numeric, mean gives the share of ones
    shareD1(i) = mean( T.d(mask) == 1 );
end

% assemble into a table
summaryTbl = table( productIDs, meanY, shareD1, ...
    'VariableNames',{'productID','meanY','shareD_d_eq_1'} )



%% 

s_omega = theta_hat(12:14); 
alpha = theta_hat(4:6); 
delta = theta_hat(7:9); 
beta = theta_hat(1:3); 
gamma0 = theta_hat(10); 
gamma1 = theta_hat(11); 
J = 3; 
 

Tsim = IE11_gen_data(T, J, s_omega, alpha, delta, beta, gamma0, gamma1); 
productIDs = unique(Tsim.productID);
nProd      = numel(productIDs);

meanY_sim   = zeros(nProd,1);
shareD1_sim = zeros(nProd,1);

for i = 1:nProd
    pid         = productIDs(i);
    mask        = Tsim.productID == pid;
    meanY_sim(i)   = mean( Tsim.Y_vj(mask) );
    shareD1_sim(i) = mean( Tsim.isPositive(mask) );
end

disp('Summary statistics of simulated data using estimated params. ')
summarySim = table( productIDs, meanY_sim, shareD1_sim, ...
    'VariableNames',{'productID','meanY','shareD_d_eq_1'} );
disp(summarySim);

