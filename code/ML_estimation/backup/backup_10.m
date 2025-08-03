%% Voyage likelihood

%%% in terms of voyage likelihood I tried with this three functions, 
% voyage lik0 is a bit more efficient than voyagelik2, and by far the slowest is voyagelik which is not optimized, 
% where I used them in the captain level likelihood and the last number of c1lik1x denotes which
% function I am using 

% Importantly they always produced the same output



% 
%  c1_lik10 =
% 
%   7.9830e-123
% 
% Elapsed time is 0.054248 seconds.
% 
% c1_lik11 =
% 
%   7.9830e-123
% 
% Elapsed time is 1.571319 seconds.
% 
% c1_lik12 =
% 
%   7.9830e-123
% 
% Elapsed time is 0.045785 seconds.
% 
% c1_lik10 =
% 
%    2.4419e-79
% 
% Elapsed time is 0.045190 seconds.
% 
% c1_lik11 =
% 
%    2.4419e-79
% 
% Elapsed time is 2.148817 seconds.
% 
% c1_lik12 =
% 
%    2.4419e-79
% 
% Elapsed time is 0.067770 seconds.
% 
% c1_lik10 =
% 
%    2.5135e-87
% 
% Elapsed time is 0.046064 seconds.
% 
% c1_lik11 =
% 
%    2.5135e-87
% 
% Elapsed time is 2.018065 seconds.
% 
% c1_lik12 =
% 
%    2.5135e-87
% 
% Elapsed time is 0.072434 seconds.
% 
% c1_lik10 =
% 
%    6.0860e-73
% 
% Elapsed time is 0.043665 seconds.
% 
% c1_lik11 =
% 
%    6.0860e-73
% 
% Elapsed time is 2.167581 seconds.
% 
% c1_lik12 =
% 
%    6.0860e-73
% 
% Elapsed time is 0.053671 seconds.
% 
% c1_lik10 =
% 
%   7.4181e-144
% 
% Elapsed time is 0.039495 seconds.
% 
% c1_lik11 =
% 
%   7.4180e-144
% 
% Elapsed time is 1.627138 seconds.
% 
% c1_lik12 =
% 
%   7.4181e-144
% 
% Elapsed time is 0.057954 seconds.

function L_v = voyageLik(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk )
    J = 3;

    % Unpack parameters
    beta        = theta(1:2);
    alpha       = theta(3:5);
    delta       = theta(6:8);
    gamma0      = theta(9);
    gamma1      = theta(10);
    sigma_omega = theta(11:13);

    L_v = 1; % initialize likelihood 

    for j = 1: J 

        if d_v(j) == 1 % positive production 
            w_hat = log(Y_v(j)) /alpha(j) - log(tau_v); 
            pr = 1 / (1 + exp(gamma0 - gamma1 * w_hat));
            f_j = normpdf(w_hat - x_v' * beta - delta(j)* a_c, 0, sigma_omega(j)); 
            L_v = L_v * pr * f_j /(alpha(j)*Y_v(j)); 

        else 
           mu    = x_v.'*beta      + delta(j)*a_c;    % scalar
           sigma = sigma_omega(j);  

           f = @(w) (1 - 1./(1 + exp(gamma0 - gamma1.*w))) .* normpdf(w, mu, sigma); 

           Pr_zero = integral(f, -Inf, Inf, 'ArrayValued', true, 'RelTol',1e-8, 'AbsTol',1e-12);
           %Pr_zero = integral(f, -Inf, Inf); 

           L_v = L_v * Pr_zero; 
        end

    end 

end

function L_v = voyageLik0(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk)

    
    % unpack
    beta        = theta(1:2);
    alpha       = theta(3:5);
    delta       = theta(6:8);
    gamma0      = theta(9);
    gamma1      = theta(10);
    sigma_omega = theta(11:13);
    
    % common term mu_j = x_v'*beta + delta_j * a_c
    mu = (x_v' * beta) + delta * a_c;      % 3×1
    
    % ---- positive‐output branch ----
    pos    = (d_v == 1);
    if any(pos)
        y_pos      = Y_v(pos);
        alpha_pos  = alpha(pos);
        sigma_pos  = sigma_omega(pos);
        mu_pos     = mu(pos);
        
        % change-of-vars
        w_hat      = log(y_pos) ./ alpha_pos - log(tau_v);       % K_pos×1  
        
        % Pr(d=1 | w)
        pr_pos     = 1 ./ (1 + exp(gamma0 - gamma1 .* w_hat));   % K_pos×1  
        
        % density in omega‐space
        f_pos      = normpdf(w_hat - mu_pos, 0, sigma_pos);      % K_pos×1  
        
        % Jacobian |dw_hat/dy|
        jac        = 1 ./ (alpha_pos .* y_pos);                  % K_pos×1  
        
        L_pos = prod(pr_pos .* f_pos .* jac);
    else
        L_pos = 1;
    end
    
    % ---- zero‐output branch ----
    zero   = ~pos;
    if any(zero)
        mu_z      = mu(zero);             % K_zero×1
        sigma_z   = sigma_omega(zero);    % K_zero×1
        
        % build a (K_zero × K) matrix of w-nodes
        %  w_jk = mu_j + sqrt(2)*sigma_j * xk_k
        W = mu_z(:) + sqrt(2) * sigma_z(:) .* xk(:).';   % size = [#zero × K]
        
        % Pr(d=1|w) at each node
        P = 1 ./ (1 + exp(gamma0 - gamma1 .* W));        % same size
        
        % ∫ [1−Pr(d=1|w)] φ((w−μ)/σ)/σ dw ≈ (1/√π)·Σ w_k·[1−P]
        Pr_zero = (1/sqrt(pi)) * ( (1 - P) * wk(:) );     % gives vector length #zero
        
        L_zero = prod(Pr_zero);
    else
        L_zero = 1;
    end
    
    % final voyage‐level likelihood
    L_v = L_pos * L_zero;
end


function L_v = voyageLik2(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk)
    % VECTORISED per-voyage likelihood using Hermite–Gauss for zeros

    % unpack parameters
    beta        = theta(1:2);
    alpha       = theta(3:5);
    delta       = theta(6:8);
    gamma0      = theta(9);
    gamma1      = theta(10);
    sigma_omega = theta(11:13);
    J = 3;

    % precompute common terms
    mu      = (x_v' * beta) + delta * a_c;          % J×1
    w_hat   = log(Y_v) ./ alpha - log(tau_v);       % J×1 (valid for Y_v>0)
    p_val   = 1 ./ (1 + exp(gamma0 - gamma1 .* w_hat));  % J×1
    f_j     = normpdf(w_hat - mu, 0, sigma_omega);  % J×1
    jac     = 1 ./ (alpha .* Y_v);                  % J×1

    L_j = zeros(J,1);
    pos_idx  = (d_v == 1);
    zero_idx = ~pos_idx;

    % positive-output branch
    if any(pos_idx)
        L_j(pos_idx) = p_val(pos_idx) .* f_j(pos_idx) .* jac(pos_idx);
    end

    % zero-output branch via Hermite–Gauss
    if any(zero_idx)
        % build matrix of nodes for each product with zero output
        sigma_z = sigma_omega(zero_idx);
        mu_z    = mu(zero_idx);
        % nodes: W_jk = mu_j + sqrt(2)*sigma_j * xk_k
        W = mu_z(:) + sqrt(2) * sigma_z(:) .* xk(:).';  % size [#zero × K]
        P = 1 ./ (1 + exp(gamma0 - gamma1 .* W));        % same size
        % integrate 1 - P over normal pdf: (1/sqrt(pi)) * Σ (1-P) * wk
        Pr_zero = (1/sqrt(pi)) * ((1 - P) * wk(:));      % #zero × 1
        L_j(zero_idx) = Pr_zero;
    end

    % voyage-level likelihood
    L_v = prod(L_j);
end


%the captain likelihood was of the form: 
function L_c = captainLik1_0(theta, d_c, Y_c, X_c, Tau_c, xk, wk)
    % CAPTAIN-LEVEL likelihood via direct integration (Eq. 8) with sigma_a = 1
    % theta(1:13) = [beta; alpha; delta; gamma0; gamma1; sigma_omega]

    % number of products and voyages
    J = 3;
    N = numel(Y_c) / J;

    % reshape data by voyage
    d_mat = reshape(d_c, J, N);        % J×N indicators
    Y_mat = reshape(Y_c, J, N);        % J×N outputs
    Tau_v = Tau_c(1:J:end);            % 1×N durations
    X_v   = X_c(1:J:end, :);           % N×p covariates

    % integrand over a_c ~ N(0,1)
    integrand = @(a_c) arrayfun(@(aa) ...
        prod(arrayfun(@(v) ...
            voyageLik0(theta, aa, d_mat(:,v), Y_mat(:,v), X_v(v,:)', Tau_v(v), xk, wk), 1:N)) ...
        .* normpdf(aa, 0, 1), ...
    a_c);

    % integrate likelihood over random effect a_c
    L_c = integral(integrand, -Inf, Inf, 'ArrayValued', true);
end

%% captain likelihood 


%%% in terms of captain likelihood I tried with  three functions, the one
%%% that has its own file and the two that are below, the
% captainlik is a bit more efficient than captainlik2, and by far the slowest is captainlik1 which is not optimized, 
% and just integrates over voyage like. 

% the one without number and 2 always produce the exact same number (they use GH integration) whereas captainlik1 uses other integration but is very close in terms 
% of results 

%some resutls: 
%Elapsed time is 0.013939 seconds.
% 
% c1_lik1 =
% 
%    1.4612e-59
% 
% Elapsed time is 0.035114 seconds.
% Elapsed time is 0.005828 seconds.
% Elapsed time is 0.001746 seconds.
% 
% c1_lik1 =
% 
%    5.2733e-92
% 
% Elapsed time is 0.037099 seconds.
% Elapsed time is 0.006170 seconds.
% Elapsed time is 0.003077 seconds.
% 
% c1_lik1 =
% 
%    2.2472e-98
% 
% Elapsed time is 0.035895 seconds.
% Elapsed time is 0.006005 seconds.
% Elapsed time is 0.006246 seconds.
% 
% c1_lik1 =
% 
%    1.0889e-67
% 
% Elapsed time is 0.063703 seconds.
% Elapsed time is 0.029792 seconds.




function L_c = captainLik1(theta, d_c, Y_c, X_c, Tau_c, xk, wk)
    % CAPTAIN-LEVEL likelihood via direct integration (Eq. 8) with sigma_a = 1
    % theta(1:13) = [beta; alpha; delta; gamma0; gamma1; sigma_omega]

    % number of products and voyages
    J = 3;
    N = numel(Y_c) / J;

    % reshape data by voyage
    d_mat = reshape(d_c, J, N);        % J×N indicators
    Y_mat = reshape(Y_c, J, N);        % J×N outputs
    Tau_v = Tau_c(1:J:end);            % 1×N durations
    X_v   = X_c(1:J:end, :);           % N×p covariates

    % integrand over a_c ~ N(0,1)
    integrand = @(a_c) arrayfun(@(aa) ...
        prod(arrayfun(@(v) ...
            voyageLik0(theta, aa, d_mat(:,v), Y_mat(:,v), X_v(v,:)', Tau_v(v), xk, wk), 1:N)) ...
        .* normpdf(aa, 0, 1), ...
    a_c);

    % integrate likelihood over random effect a_c
    L_c = integral(integrand, -Inf, Inf, 'ArrayValued', true);
end



function L_c = captainLik2(theta, d_c, Y_c, X_c, Tau_c, xk, wk)
    % CAPTAIN-LEVEL likelihood via Gauss–Hermite integration (Eq. 8)

    % scaled nodes and weights
    a_k = sqrt(2) * 1 * xk(:);       % K×1 nodes for a_c
    w_k = wk(:) / sqrt(pi);                % K×1 weights

    % determine voyages per captain
    J  = numel(theta(3:5));
    N  = numel(Y_c) / J;                    % number of voyages

    % reshape data by voyage
    d_mat   = reshape(d_c, J, N);
    Y_mat   = reshape(Y_c, J, N);
    Tau_v   = Tau_c(1:J:end);
    X_v     = X_c(1:J:end, :);

    % preallocate matrix of voyage likelihoods at each node
    L_mat = zeros(N, numel(a_k));

    % loop over Gauss–Hermite nodes
    for k = 1:numel(a_k)
        a_c = a_k(k);
        % compute each voyage's likelihood
        for v = 1:N
            d_v   = d_mat(:, v);
            Y_v   = Y_mat(:, v);
            x_v   = X_v(v, :)';
            tau_v = Tau_v(v);
            L_mat(v,k) = voyageLik0(theta, a_c, d_v, Y_v, x_v, tau_v, xk, wk);
        end
    end

    % integrate over a_c: product across voyages, then weighted sum
    prod_L = prod(L_mat, 1);                 % 1×K
    L_c    = prod_L * w_k;                   % scalar
end

