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

res = L_v_corr_v3(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);

%% 



mask  = (c_id==284);
d_cap    = d(mask);
Y_cap    = Y(mask);
Xmat_cap    = Xmat(mask,:);
tau_v_cap  = Tau(mask);

LogLc = L_c_corr_int_v3(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, ... 
                xk, wk, xk2, wk2, xk3, wk3);        

global_lik = globalLik_corr_v2(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 

%%
num_inits = 5; 
%initialize matrices to store the results 
all_theta2   = nan(numel(theta0),    2*(num_inits+1));
all_fvals2    = nan(1,        2*(num_inits+1));


theta_chol0 = to_chol_theta(theta0); 

negLL_red = @(theta_red) -globalLik_corr_v2(to_theta(theta_red), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',20);


[theta_chol_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );
theta_hat = to_theta(theta_chol_hat); % recover the full vector 

all_theta2(:,1) = theta_hat;
all_theta2(1) = fval; 



%% Loop for different initial guess. 



% Phase 1: explore relative noise scales to pick a good perturbation
noise_grid  = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1];  %% CHANGE: define relative noise grid
valid_counts = zeros(size(noise_grid));

for g = 1:numel(noise_grid)
    sigma_rel = noise_grid(g);                                               %% CHANGE: relative scale
    rel_noise = sigma_rel * randn(numel(theta_chol_hat), num_inits);         %% CHANGE: multiplicative perturbations
    inits     = repmat(theta_chol_hat,1,num_inits) .* (1 + rel_noise);       %% CHANGE: apply relative noise
    ok = false(1,num_inits);
    for j = 1:num_inits
        ok(j) = isfinite( negLL_red(inits(:,j)) );
    end
    valid_counts(g) = sum(ok);
    fprintf("σ_rel = %.3f → %2d/%2d valid inits", sigma_rel, valid_counts(g), num_inits);
end

% Choose a relative noise level with enough valid inits (e.g., first with ≥8)
good = find(valid_counts >= 8,1,'first');
if isempty(good)
    sigma_pick = noise_grid(end);
else
    sigma_pick = noise_grid(good);
end
fprintf("Using σ_rel = %.3f for perturbations", sigma_pick);




% set up
p        = numel(theta_chol0);         
theta0_mat = zeros(p, num_inits);
rng(123);                     % for reproducibility

% create 10 diverse initial guesses
for k = 1:num_inits
    % draw multiplicative factors ∼ Uniform(0.5,1.5)
    mult = 0.5 + 5*rand(p,1);
    theta0_mat(:,k) = theta_chol_hat .* mult;
end

%%
% preallocate results struct
A(num_inits) = struct('theta_hat',[],'fval',[],'exitflag',[],'output',[]);


% loop over initials
for i = 1:num_inits
    th0 = theta0_mat(:,i);
    [theta_hat, fval, exitflag, output] = fminunc(negLL_red, th0, options);
    A(i).theta_hat = theta_hat;
    A(i).fval      = fval;
    A(i).exitflag  = exitflag;
    A(i).output    = output;

    all_theta2(:, i+1) = to_theta(theta_hat); 
    all_fvals2(i+1) = fval; 

end

save('est_corr_unrestricted', 'all_theta2', 'all_fvals2');
