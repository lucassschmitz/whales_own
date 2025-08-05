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
beta   = [.0019, 0; 0.021, 0; 0.029, 0 ];        % coefficient of ship chars (weight and type).
gamma0 = 4;                % intercept for zero‐vs‐positive
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
in_somega = [1, 0.5, 0; 0.5, 3, 0; 0, 0, .5];       % initial σ_{ω,j} (σ=1)

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
% res2 = L_v_corr_v2(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res21 = exp(L_v_corr_v3(theta0, a_c, d_v, Y_v, x_v, tau_v, xk, wk, xk2, wk2, xk3, wk3) ); 
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
% res7 = L_c_corr_int_v2(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic
% res71 = L_c_corr_int_v3(theta, d_c, Y_c, x_c, tau_c, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% abs(res6- res7) < 1e-10    & abs(res6- res71) < 1e-10 
% 
% tic
% res8 = globalLik_corr(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% tic 
% [res9, B] = globalLik_corr_v2(theta0, d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3); 
% toc
% abs(res8- res9) < 1e-10


%% Estimation 

num_inits = 5; % different starting points 
negLL_red = @(chol_theta) -globalLik_corr_v2(to_theta(chol_theta), d, Y, Xmat, Tau, c_id, xk, wk, xk2, wk2, xk3, wk3 );

%matrices to store the results. 
all_theta   = nan(numel(theta0),    num_inits+1);
all_fvals    = nan(1,            num_inits+1);

% % Diagnostic1
% % Add this after creating theta_chol0
% fprintf('--- Starting Diagnostic ---\n');
% initial_negLL = negLL_red(theta_chol0);
% 
% fprintf('Initial objective function value: %f\n', initial_negLL);
% if isnan(initial_negLL) || isinf(initial_negLL)
%     fprintf('WARNING: Objective function returned NaN or Inf!\n');
% end
% fprintf('--- End Diagnostic ---\n');
% 
% theta_chol02 = theta_chol0; 
% theta_chol02(1:3) = [.2; .2; .2];
% initial_negLL2 = negLL_red(theta_chol02);
%
% %chekc that objective is not flat 
% a = negLL_red(theta_chol0) 
% b = negLL_red(theta_chol0+1e-4)


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
    'MaxIterations',200, ... 
    'MaxFunctionEvaluations', 20, ...
    'OptimalityTolerance',1e-20,...
    'StepTolerance',      1e-15) 
 
[theta_red_hat, fval, exitflag, output] = fminunc(negLL_red, theta_chol0, options );


theta_hat = to_theta(theta_red_hat); % recover the full vector 
all_theta(:,1) = theta_hat;
all_theta(1) = fval; 
 

%% Initial values
 

% Store results in a table for easier comparison
results_table = table('Size', [num_inits, 4], 'VariableTypes', {'double', 'cell', 'cell', 'double'}, ...
                      'VariableNames', {'InitialValueID', 'EstimatedTheta', 'DifferenceFromReal', 'Distance'});

% Generate different initial values by adding noise directly to the Cholesky parameters
theta_chol_inits = repmat(theta_red_hat, 1, num_inits) + randn(size(theta_red_hat, 1), num_inits) * 0.01;


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

save('est_corr_unrestricted', 'all_theta', 'all_fvals');




%% Comparison of Estimates with theta_real

fprintf('\n\n--- Comparison of Final Estimates with True Theta ---\n\n');
fprintf('Theta Real:\n');
disp(theta_real');
disp('--- Estimation Results ---');
disp(results_table);

