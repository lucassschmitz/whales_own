clear
clc

rng(123)
% parameters 
C = 55;       % Number of captains 
J = 3;         % 3 product‐types: bones, oil, sperm
Vmax = 30;     % maximum of captain voyages. 

s_omega = 2 * eye(J);     % std. dev. of product‐specific noise omega_vj
s_omega = [2 1 1; 1 2 1; 1 1 2];
disp(['Lowest eigenvalue = ' num2str(min(eig(s_omega)))]);

 
alpha  = [1.3; .8; .82];    % exponent 
delta  = [2; 2.3; 2.9];    % coefficient of captain random effects
beta   = [.019, 0.3; 0.021, 0; 0.029, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 4;                % intercept for zero‐vs‐positive
gamma1 = 2;                % slope on log(w_{1vj})


Tsim = simulate_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1);

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
in_alpha     = [.8;.8; .8];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %6:8
in_gamma0 = .4;               % initial γ0 9 
in_gamma1 = 1;               % initial γ1 10
in_somega = [1, 0.5, 0; 0.5, 3, 0; 0, 0, .5];       % initial σ_{ω,j} (σ=1)

% parameter guess temp 
in_beta   = [0.1; 0.1; 0.1];          % initial β's
in_alpha     = [.8;.8; .8];       % initial   α's (α=1)
in_delta = [1 ; 1; 1]; %6:8
in_gamma0 = .4;               % initial γ0 9 
in_gamma1 = 1;               % initial γ1 10
in_somega = [10,0, 0; 0, 30, 0; 0, 0, 10];       % initial σ_{ω,j} (σ=1)

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:)];
theta_chol0 = to_chol_theta(theta0); % theta reduced 

[xk, wk] = HG(15);  % 30-node rule
[xk2, wk2] = HG2D(15); 
[xk3, wk3] = HG3D(15);

% CODE USED TO CHECK THE FUNCTIONS 
% a_c = .76; d_v = d(1:3); Y_v = Y(1:3); x_v = Xmat(3,:)'; tau_v = Tau(3); 
% a_c = .76; d_v = d(4:6); Y_v = Y(4:6); x_v = Xmat(4,:)'; tau_v = Tau(4); 
% a_c = .76; d_v = d(19:21); Y_v = Y(19:21); x_v = Xmat(19,:)'; tau_v = Tau(19); 
% a_c = .76; d_v = d(28:30); Y_v = Y(28:30); x_v = Xmat(28,:)'; tau_v = Tau(28); 
% a_c = .76; d_v = d(46:48); Y_v = Y(46:48); x_v = Xmat(46,:)'; tau_v = Tau(28); 
% a_c = .76; d_v = d(157:159); Y_v = Y(157:159); x_v = Xmat(157,:)'; tau_v = Tau(157); 
%
% the 2nd and 3rd versions of L_v_corr can be order of magnitude faster
% than the 1st version, spetially if the voyage has zero production for the
% three porducts, in which case the computing times can be 48s first
% version and , .00004s for the other ones. 
% tic
% res = L_v_corr(theta0, a_c, d_v, Y_v, x_v, tau_v); 
% toc
% tic
% res2 = L_v_corr(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res21 = exp(L_v_corr(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3) ); 
% toc
% theta2 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; diag(in_somega)];
% res3 = voyageLik(theta2, a_c, d_v, Y_v, x_v, tau_v, xk, wk); 
%
% mask   = (c_id == 3);
% d_c    = d(mask);
% Y_c    = Y(mask);
% x_c    = Xmat(mask, :);
% tau_c  = Tau(mask);
% a_c = 0; 
% theta = theta0; 
% tic
% res5 = L_c_corr(theta0, a_c, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res6 = L_c_corr_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res7 = L_c_corr_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res71 = L_c_corr_int(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% abs(res6- res7) < 1e-10    & abs(res6- res71) < 1e-10 
% 
% tic
% res8 = globalLik_corr(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic 
% [res9, B] = globalLik_corr(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% abs(res8- res9) < 1e-10


%% Estimation 

num_inits = 2; % different starting points 
negLL_red = @(chol_theta) -globalLik_corr(to_theta(chol_theta), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

%matrices to store the results. 
all_theta   = nan(numel(theta0),    num_inits+1);
all_fvals    = nan(1,            num_inits+1);

% % Diagnostic1
initial_negLL = negLL_red(theta_chol0);

fprintf('Initial objective function value: %f\n', initial_negLL);
if isnan(initial_negLL) || isinf(initial_negLL)
    fprintf('WARNING: Objective function returned NaN or Inf!\n');
end
fprintf('--- End Diagnostic ---\n');

theta_chol02 = theta_chol0; 
theta_chol02(1:3) = [.2; .2; .2];
initial_negLL2 = negLL_red(theta_chol02);

%chekc that objective is not flat 
a = negLL_red(theta_chol0) 
b = negLL_red(theta_chol0+1e-4)


% % --- Diagnostic2: Plot Likelihood Surface ---
% param_to_test = 1; % Test the first parameter
% theta_start = theta_chol0;
% n_points = 50;
% param_range = linspace(theta_start(param_to_test) * 0.9, theta_start(param_to_test) * 1.1, n_points);
% y_vals = zeros(n_points, 1);
% 
% for i = 1:n_points
%     theta_test = theta_start;
%     theta_test(param_to_test) = param_range(i);
%     y_vals(i) = negLL_red(theta_test);
% end
% 
% figure;
% plot(param_range, y_vals);
% title(['Log-Likelihood vs. Parameter ', num2str(param_to_test)]);
% xlabel(['Value of Parameter ', num2str(param_to_test)]);
% ylabel('Negative Log-Likelihood');

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',300, ... 
    'MaxFunctionEvaluations', 20, ...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',      1e-15) 
 
[theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );


theta_hat = to_theta(theta_red_hat); % recover the full vector 
all_theta(:,1) = theta_hat;
all_fvals(1) = fval; 
 

%% Initial values
 

% Store results in a table for easier comparison
results_table = table('Size', [num_inits, 4], 'VariableTypes', {'double', 'cell', 'cell', 'double'}, ...
                      'VariableNames', {'InitialValueID', 'EstimatedTheta', 'DifferenceFromReal', 'Distance'});

% 1. choose how many inits and what noise‐levels to try
noise_grid  = [0.1, 0.02, 0.1, 0.2];
valid_counts = zeros(size(noise_grid)); % much noise generates NaN 

for g = 1:length(noise_grid)
    sigma = noise_grid(g);
    % generate your inits for this σ
    rel_noise = sigma * randn(numel(theta_red_hat), num_inits);
    inits     = repmat(theta_red_hat,1,num_inits) .* (1 + rel_noise);

    % test each one
    ok = false(1, num_inits);
    for j = 1:num_inits
        val = negLL_red(inits(:,j));
        ok(j) = isfinite(val);
    end

    valid_counts(g) = sum(ok);
    fprintf('σ = %.3f  →  %2d / %2d valid inits\n', sigma, valid_counts(g), num_inits);
end


% 2) Once you pick σ_rel, rebuild inits and run your minimizations:
sigma_rel       = 0.9;  % for example
rel_noise   = sigma_rel * randn(numel(theta_red_hat), num_inits);
theta_chol_inits = repmat(theta_red_hat,1,num_inits) .* (1 + rel_noise);




%% 
for i = 1:num_inits
    fprintf('Running Minimization for Initial Value %d...\n', i); 

    % if aux == NaN does not use the intiial value to minimize.  
    aux = negLL_red(theta_chol_inits(:,i))
    if isnan(aux) || isinf(aux)
        fprintf('Initial value %d returned NaN or Inf. Skipping.\n', i);
        continue;
    end    

    % Run the optimizer
    [theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol_inits(:, i), options );
    
    % Store results
    theta_hat = to_theta(theta_red_hat);
    difference = theta_hat - theta_real;

    all_theta(:,   i+1) = theta_hat;
    all_fvals(  i+1)  = fval;

    distance = norm(difference);
    
    results_table.InitialValueID(i) = i;
    results_table.EstimatedTheta{i} = theta_hat';
    results_table.DifferenceFromReal{i} = difference';
    results_table.Distance(i) = distance;
end

save('output/sim_corr_unrestricted', 'all_theta', 'all_fvals');




%% Comparison of Estimates with theta_real

fprintf('\n\n--- Comparison of Final Estimates with True Theta ---\n\n');
fprintf('Theta Real:\n');
disp(theta_real');
disp('--- Estimation Results ---');
disp(results_table);

%% ==== Post-Estimation Analysis for sim_correlation ====

% Load the estimation results from the simulation
load('output/sim_corr_unrestricted.mat', 'all_theta', 'all_fvals');

% --- 1) Filter out invalid runs (NaN fvals) ---
validCols = ~isnan(all_fvals);
theta_mat = all_theta(:, validCols);
fvals = all_fvals(validCols);

% --- 2) Summary statistics for each parameter ---
mean_theta = mean(theta_mat, 2)';
std_theta = std(theta_mat, 0, 2)';
var_theta = var(theta_mat, 0, 2)';
cov_theta = cov(theta_mat');
cv_theta = std_theta ./ abs(mean_theta);

paramNames = {'beta1','beta2','beta3', 'alpha1','alpha2','alpha3', ...
              'delta1','delta2','delta3', 'gamma0','gamma1', ...
              'sigma11','sigma12','sigma13','sigma21','sigma22', ...
              'sigma23','sigma31','sigma32','sigma33'}';

summaryStats = table(mean_theta', std_theta', var_theta', cv_theta', ...
    'VariableNames', {'Mean','StdDev','Variance','CV'}, 'RowNames', paramNames);

disp('Parameter summary statistics from simulation (valid runs only):');
disp(summaryStats);

% --- 3) Choose the best theta and compare with the true parameters ---
[~, relIdx] = min(fvals);
validIdxs = find(validCols);
bestIdx = validIdxs(relIdx);
theta_best = all_theta(:, bestIdx);

% Comparison with the known "real" theta from the simulation
bestTable = table(theta_real, theta_best, (theta_best - theta_real), ...
    'VariableNames', {'Theta_Real', 'Theta_Best_Estimated', 'Difference'}, ...
    'RowNames', paramNames);

disp('Best estimated parameters vs. true parameters:');
disp(bestTable);

% --- 4) Simulate data with theta_best and compare moments ---
% Unpack the best-fitting parameters
beta_est = theta_best(1:3);
alpha_est = theta_best(4:6);
delta_est = theta_best(7:9);
gamma0_est = theta_best(10);
gamma1_est = theta_best(11);
Sigma_omega_est = reshape(theta_best(12:end), 3, 3);

% --- "Real" moments from the initial simulated data (Tsim) ---
productIDs = unique(Tsim.productID);
real_meanY = zeros(J,1);
real_share = zeros(J,1);
real_med = zeros(J,1);

for i = 1:J
    pid = productIDs(i);
    mask = Tsim.productID == pid;
    product_data = Tsim.Y_vj(mask);
    pos_product = product_data(product_data > 0);
    
    real_meanY(i) = mean(product_data);
    real_share(i) = mean(Tsim.isPositive(mask));
    real_med(i) = median(pos_product);
end

% --- Simulate new data using the best estimated parameters ---
nSims = 40;
sim_meanY = zeros(J, nSims);
sim_medianY = zeros(J, nSims);
sim_share = zeros(J, nSims);
simY_new = cell(J,1);
for j = 1:J
    simY_new{j} = [];
end


% --- FIX: Rename Nvoy to voyageID to match the function's expectation ---
Tsim.Properties.VariableNames{'Nvoy'} = 'n_voy';
Tsim.Properties.VariableNames{'Duration'} = 'Tau';
[Tsim.voyageID, ~] = findgroups(Tsim.captainID, Tsim.n_voy);

for s = 1:nSims
    % Use IE11_gen_data to keep covariates (tonnage, duration) from Tsim constant
    Tsim_new = IE11_gen_data(Tsim, J, Sigma_omega_est, alpha_est, delta_est, beta_est, gamma0_est, gamma1_est);
    for i = 1:J
        pid = productIDs(i);
        m = Tsim_new.productID == pid;
        prod_new = Tsim_new.Y_vj(m);
        pos_prod_new = prod_new(prod_new > 0 & isfinite(prod_new));
        
        sim_meanY(i,s) = mean(prod_new);
        sim_medianY(i,s) = median(pos_prod_new);
        sim_share(i,s) = mean(Tsim_new.isPositive(m));
        simY_new{i} = [simY_new{i}; prod_new];
    end
end

% --- Average simulated moments ---
meanY_sim = mean(sim_meanY, 2);
medY_sim = mean(sim_medianY, 2);
share_sim = mean(sim_share, 2);

% --- Comparison Table ---
compareTbl = table(productIDs, real_meanY, meanY_sim, real_share, share_sim, real_med, medY_sim, ...
    'VariableNames', {'ProductID','Real_MeanY','Sim_MeanY','Real_Share', ...
                    'Sim_Share','Real_Median','Sim_Median'});
disp('"Real" vs. Simulated moments (median conditional on positive production):');
disp(compareTbl);

% --- 5) Histograms of "real" vs. newly simulated output ---
figure;
for j = 1:J
    pid = productIDs(j);
    
    % Prepare "Real" (original Tsim) data
    Yr = Tsim.Y_vj(Tsim.productID == pid);
    Yr_pos = Yr(Yr > 0 & isfinite(Yr));
    
    % Prepare newly simulated data
    Ys = simY_new{j};
    Ys_pos = Ys(Ys > 0 & isfinite(Ys));

    % Plot "Real" Data
    subplot(2, J, j);
    histogram(Yr_pos, 'Normalization', 'pdf');
    xlim([0, prctile(Yr_pos, 99)]);
    title(sprintf('Original Sim: Product %d', pid));
    xlabel('Y');
    ylabel('Density');

    % Plot Newly Simulated Data
    subplot(2, J, J + j);
    histogram(Ys_pos, 'Normalization', 'pdf');
    xlim([0, prctile(Ys_pos, 99)]);
    title(sprintf('New Sim (from est.): Product %d', pid));
    xlabel('Y');
    ylabel('Density');
end