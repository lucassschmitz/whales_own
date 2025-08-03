%% minimization stability 
 
clc 
clear
%── 1) Load saved estimates ──────────────────────────────────────────────
load('est_corr_unrestricted.mat','theta_hat')    % fminunc initial
theta_hat = from_theta(theta_hat); 


%load('est_corrind__restricted.mat','theta_hat2')  % fmincon initial
load('est_corr_unrestrictedMat.mat','A')    % struct array from fminunc runs
%load('est_ind__restrictedMat.mat','A2')  % struct array from fmincon runs

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


%%

% fmincon block
% base = N+1;
% allTheta(base+1,   :) = theta_hat2';
% for i = 1:N
%     allTheta(base+1+i, :) = A2(i).theta_hat';
% end

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
    'sigma_bone','sigma_sperm','sigma_oil', 's_12', 's_13', 's_23' ...    % σ_ω parameters
}'; %paramNames = compose("theta%02d", 1:p)';  
params_tab = table(mean_theta', std_theta', var_theta', cv_theta', ...
          'VariableNames', {'Mean','StdDev','Variance','CV'}, ...
          'RowNames', paramNames);
disp(params_tab)





 function theta_red = from_theta(theta)
    % FROM_THETA  Extracts the 17 unique parameters from full 20×1 theta
    %  theta: [β(3); α(3); δ(3); γ0; γ1; vec(Sigma_omega) (9)]
    Sigma = reshape(theta(12:20), 3, 3);  % full 3×3 covariance matrix
    theta_red = [
        theta(1:11);      % β1–3, α1–3, δ1–3, γ0, γ1
        Sigma(1,1);       % σ11
        Sigma(2,2);       % σ22
        Sigma(3,3);       % σ33
        Sigma(1,2);       % σ12 = σ21
        Sigma(1,3);       % σ13 = σ31
        Sigma(2,3)        % σ23 = σ32
    ];
 end
