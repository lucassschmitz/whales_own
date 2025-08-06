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
in_lambda = [5;  1; 1]; 
% σ_a is normalized to 1 

theta0 = [in_beta; in_alpha; in_delta; in_gamma0; in_gamma1; in_somega(:); in_lambda];

[xk, wk] = HG(30);  % 30-node rule
[xk2, wk2] = HG2D(15); 
[xk3, wk3] = HG3D(15);

a_c = .76; 
d_v = d(1:3); 
Y_v = Y(1:3); 
x_v = Xmat(3,:)';
tau_v = Tau(3); 

res = L_v_corr_v4(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3);


mask  = (c_id==284);
d_cap    = d(mask);
Y_cap    = Y(mask);
Xmat_cap    = Xmat(mask,:);
tau_v_cap  = Tau(mask);

LogLc = L_c_corr_int_v4(theta0, d_cap, Y_cap, Xmat_cap, tau_v_cap, ... 
                xk, wk, xk2, wk2, xk3, wk3);        


global_lik = globalLik_corr_v3(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 

%%
num_inits = 3; 
%initialize matrices to store the results 
all_theta2   = nan(numel(theta0),    2*(num_inits+1));
all_fvals2    = nan(1,        2*(num_inits+1));


theta_chol0 = to_chol_theta(theta0); 

negLL_red = @(theta_red) -globalLik_corr_v3(to_theta(theta_red), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

% Optimization options
options = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','iter', ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',10);


[theta_chol_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );
theta_hat = to_theta(theta_chol_hat); % recover the full vector 

all_theta2(:,1) = theta_hat;
all_fvals2(1) = fval; 

%% Generate alternative initial guesses 
% to add perturbations big enough we see how big the perturbations have to be to generate prroblems evaluating htem 


% Phase 1: explore relative noise scales to pick a good perturbation
noise_grid  = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1];  %% CHANGE: define relative noise grid
valid_counts = zeros(size(noise_grid));

for g = 1:numel(noise_grid)
    sigma_rel = noise_grid(g);                                               
    rel_noise = sigma_rel * randn(numel(theta_chol_hat), num_inits);       
    inits     = repmat(theta_chol_hat,1,num_inits) .* (1 + rel_noise);  
    ok = false(1,num_inits);
    for j = 1:num_inits
        ok(j) = isfinite( negLL_red(inits(:,j)) );
    end
    valid_counts(g) = sum(ok);
    fprintf("σ_rel = %.3f → %2d/%2d valid inits \n", sigma_rel, valid_counts(g), num_inits);
end

sigma_pick = .5; % can be changed. 

% Phase 2: generate perturbations and run minimizations
rel_noise   = sigma_pick * randn(numel(theta_chol_hat), num_inits);       
theta0_mat  = repmat(theta_chol_hat,1,num_inits) .* (1 + rel_noise); 


%%
% preallocate results struct
A(num_inits) = struct('theta_hat',[],'fval',[],'exitflag',[],'output',[]);


% loop over initials
for i = 1:num_inits
    th0 = theta0_mat(:,i);

    if ~isfinite(negLL_red(th0))
    fprintf('Skipping init %d (invalid)\n', i);
    continue;
    end
    [theta_hat, fval, exitflag, output] = fminunc(negLL_red, th0, options);
    A(i).theta_hat = theta_hat;
    A(i).fval      = fval;
    A(i).exitflag  = exitflag;
    A(i).output    = output;

    all_theta2(:, i+1) = to_theta(theta_hat); 
    all_fvals2(i+1) = fval; 

end

save('est_corr_unrestricted_v2', 'all_theta2', 'all_fvals2');


%% ==== Post-Estimation Analysis ====
% Load saved results

load('est_corr_unrestricted_v2.mat','all_theta2', 'all_fvals2') 


% 1) Filter out invalid runs (NaN fvals)
validCols   = ~isnan(all_fvals2);                         
theta_mat   = all_theta2(:, validCols);
fvals       = all_fvals2(validCols);

% 2) Summary statistics for each parameter (over valid runs)
mean_theta = mean(theta_mat, 2)';                        
std_theta  = std(theta_mat, 0, 2)';                       
var_theta  = var(theta_mat, 0, 2)';                       
cov_theta  = cov(theta_mat');                             
cv_theta   = std_theta ./ abs(mean_theta);                

paramNames = { 'beta1','beta2','beta3', ...
               'alpha1','alpha2','alpha3', ...
               'delta1','delta2','delta3', ...
               'gamma0','gamma1', ...
               'sigma11','sigma12','sigma13','sigma21','sigma22','sigma23','sigma31','sigma32','sigma33', 'lambda1', 'lambda2', 'lambda3'  }';

summaryStats = table(mean_theta', std_theta', var_theta', cv_theta', ...
    'VariableNames', {'Mean','StdDev','Variance','CV'}, 'RowNames', paramNames); 

disp('Parameter summary statistics (valid runs only):');
disp(summaryStats);

% 3) Choose θ_best that minimizes the negative log-likelihood among valid runs
[~, relIdx] = min(fvals);                               
validIdxs   = find(validCols);
bestIdx     = validIdxs(relIdx);                        
theta_best  = all_theta2(:, bestIdx);                   

bestTable = table(theta_best, 'VariableNames', {'theta_best'}, 'RowNames', paramNames);   

disp('Best parameter estimates:' );
disp(bestTable);

% 3) Simulate data with θ_best and compare moments to real data
% Unpack θ_best
beta   = theta_best(1:3);
alpha  = theta_best(4:6);
delta  = theta_best(7:9);
gamma0 = theta_best(10);
gamma1 = theta_best(11);
Sigma_omega = reshape(theta_best(12:20), 3, 3);
lambda = theta_best(end-2: end );


% Real data moments
productIDs = unique(T.productID);
real_meanY = zeros(J,1);
real_share = zeros(J,1);
real_med = zeros(J,1);
for i = 1:3
    pid = productIDs(i);
    mask = T.productID == pid;
    real_meanY(i) = mean(T.Y(mask));
    real_share(i)= mean(T.d(mask)==1);

    product = T.Y(mask); 
    pos_product = product(product>0); 

    real_med(i) = median(pos_product); 

end


% Simulate
nSims      = 40;
sim_meanY  = zeros(J, nSims);
sim_share  = zeros(J, nSims);
simY = cell(J,1);                            %% CHANGE
for j = 1:J                                  %% CHANGE
    simY{j} = [];
end
for s = 1:nSims
    Tsim = IE12_gen_data(T, 3, Sigma_omega, alpha, delta, beta, gamma0, gamma1, lambda);  
    for i = 1:J
        pid  = productIDs(i);
        m    = Tsim.productID == pid;
        prod = Tsim.Y_vj(m); 
        pos_prod = prod(prod>0); 
        sim_meanY(i,s) = mean(Tsim.Y_vj(m));
     
        sim_medianY(i,s) = median(pos_prod); 

        sim_share(i,s)= mean(Tsim.isPositive(m));
        simY{i} = [simY{i}; Tsim.Y_vj(m)];
    end
end

% Average simulated moments
meanY_sim   = mean(sim_meanY, 2);
medY_sim   = mean(sim_medianY, 2);


share_sim   = mean(sim_share, 2);

% Comparison table
compareTbl = table(productIDs, real_meanY, real_share,  real_med, meanY_sim, share_sim, medY_sim,   ...
    'VariableNames', {'productID', 'RealMeanY','RealShare', 'RealMedian','SimMeanY', 'SimShare', 'SimMedian'}); 
disp('Real vs Simulated moments (median conditional on positive production):');
disp(compareTbl);


%%
% 6) Histograms of real vs simulated output by product
figure;                                          
for j = 1:J                                  
    pid = productIDs(j);                         
    Yr  = T.Y(T.productID==pid);
    Yr = Yr(Yr> 0); 

    Ys  = simY{j};                      
    Ys = Ys(Ys>0); 

    %Ys = Ys(Ys< prctile(Ys, 80)); 
    
    subplot(2, J, j);                      %% CHANGE
    histogram(Yr, 'Normalization','pdf');   
    xlim([0, prctile(Yr, 99)]);
    title(sprintf('Real: Product %d positive only', pid));    %% CHANGE
    xlabel('Y'); ylabel('Density');               %% CHANGE

    subplot(2, J, J+j);                 %% CHANGE
    %histogram(Ys ,'NumBins', 50);  
    histogram(Ys, 'Normalization', 'pdf');
    xlim([0, prctile(Ys, 99)]);
    title(sprintf('Sim: Product %d positive only', pid));     %% CHANGE
    xlabel('Y'); ylabel('Density');               
end                                             %% CHANGE
 