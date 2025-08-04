%% Code to estimate the parameters of the model using the real data. 

clear
clc
rng(123)

T = readtable('../../Data/temps/clean_ML_estimation.xlsx');
disp(head(T));


%parameters 
C = numel( unique(T.captainID) );     % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = max(T.n_voy);  % maximum of captain voyages. 

%prepare data vectors
n     = size(T,1);
d     = T.d;
Y     = T.Y;
j_i   = T.productID;
Xmat  = T.X1;
Tau = T.Tau;
c_id = T.captainID; % use built‑in captain ID

%parameter guess. 
in_beta   = [0.01; 0.01; 0.01];          % initial β's
in_alpha     = [1; 1; 1];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %7:9
in_gamma0 = .4;               % initial γ0 10 
in_gamma1 = 1;               % initial γ1 11
in_somega = [1; 3; .5];       %positions 12:14
%σ_a is normalized to 1 

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];

[xk, wk] = HG(30);  % 30-node rule

a_c = .76; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 

res = voyageLik_v2(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 

mask  = (c_id==284);
d_cap    = d(mask);
Y_cap    = Y(mask);
Xmat_cap    = Xmat(mask,:);
tau_v_cap  = Tau(mask);

LogLc = captainLik_v2(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, xk, wk);        

L_global3 = globalLik_v3(theta0, d, Y, Xmat, Tau, c_id, xk, wk);

%Set up negative log-likelihood for minimization
negLL = @(th) -globalLik_v3(th, d, Y, Xmat, Tau, c_id, xk, wk);

%Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',20, ...
    'MaxFunctionEvaluations',20);

%%

% estimate theta by minimizing negative log-likelihood
[theta_hat, fval, exitflag, output] = fminunc(negLL, theta0, options);

%Final log-likelihood at estimated parameters
ll_hat = globalLik_v3(theta_hat, d, Y, Xmat, Tau, c_id, xk, wk);

save('est_ind_unrestricted', 'theta_hat')

%% Loop for different initial guess. 

% set up
nInits   = 10;                     % how many initials you want
p        = numel(theta0);         
theta0_mat = zeros(p, nInits);
rng(123);                     % for reproducibility

% create 10 diverse initial guesses
for k = 1:nInits
    % draw multiplicative factors ∼ Uniform(0.5,1.5)
    mult = 0.5 + 5*rand(p,1);
    theta0_mat(:,k) = theta0 .* mult;
end

% preallocate results struct
A(nInits) = struct('theta_hat',[],'fval',[],'exitflag',[],'output',[]);


% loop over initials
for i = 1:nInits
    th0 = theta0_mat(:,i);
    [theta_hat, fval, exitflag, output] = fminunc(negLL, th0, options);
    A(i).theta_hat = theta_hat;
    A(i).fval      = fval;
    A(i).exitflag  = exitflag;
    A(i).output    = output;
end

save('est_ind_unrestrictedMat.mat','A');


%%
theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega];

%Lower bounds for non-negativity on alpha (4:6),  delta (7:9), gamma1(11), sigma_omega(12:14)
lb = -inf(size(theta0));
lb([4:6, 7:9,11:14]) = 0;

options = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'TolFun',1e-4, 'TolX',1e-4, 'OptimalityTolerance',1e-4, ...
    'MaxIterations',10, 'MaxFunctionEvaluations',3000);

[theta_hat2, fval, exitflag, output] = fmincon(negLL, theta0, [], [], [], [], lb, [], [], options);
save('est_ind_restricted', 'theta_hat2')

A2(nInits) = struct('theta_hat',[],'fval',[],'exitflag',[],'output',[]); % preallocate struct array

for i = 1:nInits
    theta0 = theta0_mat(:, i);

    [theta_hat2, fval, exitflag, output] = fmincon(negLL, theta0, [], [], [], [], lb, [], [], options);

    %Store results
    A2(i).theta_hat = theta_hat2;   % column vector
    A2(i).fval      = fval;        % scalar
    A2(i).exitflag  = exitflag;    % scalar
    A2(i).output    = output;      % struct (e.g., iterations, message, etc.)
end

save('est_ind_restrictedMat.mat', 'A2');


%% minimization stability 
 

% 1) Load saved estimates ──────────────────────────────────────────────
load('est_ind_unrestricted.mat','theta_hat')    % fminunc initial
load('est_ind_restricted.mat','theta_hat2')  % fmincon initial
load('est_ind_unrestrictedMat.mat','A')    % struct array from fminunc runs
load('est_ind_restrictedMat.mat','A2')  % struct array from fmincon runs

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

% ── 3) Compute dispersion statistics ─────────────────────────────────────
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

%% select minimum of all the local minima
fval_unc    = negLL(theta_hat);
fval_con    = negLL(theta_hat2);
fvals_unc   = [A.fval];
fvals_con   = [A2.fval];

allFvals    = [fval_unc, fval_con, fvals_unc, fvals_con];
allThetas   = [theta_hat, theta_hat2, reshape([A.theta_hat], [], numel(A)), reshape([A2.theta_hat], [], numel(A2))];


% Find the index of the minimal negLL
[~, bestIdx] = min(allFvals);

% Select best parameter vector
bestTheta   = allThetas(:, bestIdx);


% Present best parameter estimates
best_params_tbl = table(bestTheta, ...
    'RowNames', paramNames, ...
    'VariableNames', {'bestTheta'});
disp('Best parameter estimates:');
disp(best_params_tbl);

%% model fit. 


% moments in the data 
% assume T is your table
productIDs = unique(T.productID);


nProd      = numel(productIDs);

%pre‐allocate
meanY   = zeros(nProd,1);
shareD1 = zeros(nProd,1);

for i = 1:nProd
    pid      = productIDs(i);
    mask     = T.productID == pid;
    meanY(i) = mean( T.Y(mask) );
    %if d is already 0/1 numeric, mean gives the share of ones
    shareD1(i) = mean( T.d(mask) == 1 );
end

% assemble into a table
summaryTbl = table( productIDs, meanY, shareD1, ...
    'VariableNames',{'productID','meanY','shareD_d_eq_1'} )



%%

s_omega = bestTheta(12:14); 
alpha = bestTheta(4:6); 
delta = bestTheta(7:9); 
beta = bestTheta(1:3); 
gamma0 = bestTheta(10); 
gamma1 = bestTheta(11); 
J = 3; 
 
nSims = 100; 
productIDs = unique(T.productID);
nProd      = numel(productIDs);


% Run simulations and compute moments
for s = 1:nSims
    Tsim = IE11_gen_data(T, J, s_omega, alpha, delta, beta, gamma0, gamma1);
    for i = 1:nProd
        pid  = productIDs(i);
        mask = Tsim.productID == pid;
        meanY_all(i,s)   = mean(Tsim.Y_vj(mask));
        shareD1_all(i,s) = mean(Tsim.isPositive(mask));
    end
end

% Compute average moments across all simulations
meanY_bar   = mean(meanY_all, 2);
shareD1_bar = mean(shareD1_all, 2);

% Display average summary statistics
disp('Average summary statistics across 100 simulations:');
summaryAvg = table(productIDs, meanY_bar, shareD1_bar, ...
    'VariableNames', {'productID','meanY_avg','shareD_d_eq_1_avg'});
disp(summaryAvg);


%%
alpha_alt = [1.13; 0.52; 0.63];  % new alphas
theta_alt = bestTheta; 
theta_alt(4:6) = alpha_alt; 
negLL(theta_alt) %worse than the min of allFvals. 



nSims_alt = 100; 
meanY_alt   = zeros(nProd, nSims_alt);
shareD1_alt = zeros(nProd, nSims_alt);
for s = 1:nSims_alt
    Tsim = IE11_gen_data(T, J, s_omega, alpha_alt, delta, beta, gamma0, gamma1);
    for i = 1:nProd
        pid  = productIDs(i);
        mask = Tsim.productID == pid;
        meanY_alt(i,s)   = mean(Tsim.Y_vj(mask));
        shareD1_alt(i,s) = mean(Tsim.isPositive(mask));
    end
end
meanY_bar_alt   = mean(meanY_alt, 2);
shareD1_bar_alt = mean(shareD1_alt, 2);

disp('Average moments with alternative alphas (100 sims):');
summaryAlt = table(productIDs, meanY_bar_alt, shareD1_bar_alt, ...
    'VariableNames',{'productID','meanY_avg','shareD_d_eq_1_avg'});
disp(summaryAlt);

