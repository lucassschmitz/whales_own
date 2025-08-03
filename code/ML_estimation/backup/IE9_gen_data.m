function Tsim = IE9_gen_data(C, J, Vmax, s_omega, alpha, delta, beta, gamma0, gamma1)

% the difference with IE8_gen_data is that we normalize \sigma_a = 1
% because we can always readjust the \delta to account for its scale. 

% Pre‐allocate the “long” vectors exactly as before
N  = C * Vmax * J;
captainID = zeros(N,1);
Nvoy = zeros(N,1);
productID = zeros(N,1);

X1        = zeros(N,1);
X2        = zeros(N,1);

a_c       = zeros(N,1);
omega_vj  = zeros(N,1);   % STILL a scalar per row, but filled from a J‐vector
log_wvj   = zeros(N,1);
P_pos     = zeros(N,1);
isPositive= zeros(N,1);
Y_vj      = zeros(N,1);
Dur       = zeros(N,1); 

idx = 0;




% Pre‐draw all a_c for each captain
true_ac = randn(C,1);

for c = 1:C
    for v = 1:Vmax
        

        % Draw a J‐vector of omegas for this (captain c, voyage v):
        Omega_vec = mvnrnd( zeros(1,J), s_omega );  
        
        % voyage duration 
        Dur_aux = randi([1,10]); 

        idx_start = J * (Vmax * (c-1) + (v-1)) + 1;
        idx_end = idx_start + J - 1;

        Dur(idx_start:idx_end) = Dur_aux;

        % ship chars. 
        x1_v = 10 + 0.1*randn;
        x2_v = rand > 0.5;

        % 2) Loop over products j = 1..J, assign each component of Omega_vec:
        for j = 1:J
            idx = idx + 1;

            captainID(idx) = c;
            Nvoy(idx)  = v;
            productID(idx)  = j;

            % % Covariates X1, X2
            X1(idx) = x1_v;
            X2(idx) = x2_v;

            % Captain skill, replicated for every row with the same captain c:
            a_c(idx) = true_ac(c);

            % Now assign the j‐th component of the MVN draw:
            omega_vj(idx) = Omega_vec(j);

            % Compute log w_{vj} = X_f * beta + delta_j * a_c + omega_vj:
            log_wvj(idx) = beta(1)*x1_v + beta(2)*x2_v + delta(j)*a_c(idx) + omega_vj(idx);
                          
            P_pos(idx) = 1 / (1 + exp(gamma0 - gamma1 * log_wvj(idx)));
            
            isPositive(idx) = (rand < P_pos(idx));

            tau_v = 1;  % keep duration fixed for now
            if isPositive(idx)
                Y_vj(idx) = (tau_v * exp(log_wvj(idx)))^alpha(j);
            else
                Y_vj(idx) = 0;
            end
        end
    end
end

 
% Finally, pack into a table as before:
Tsim = table( ...
   captainID, Nvoy, productID, ...
   X1, X2, ...
   a_c, omega_vj, log_wvj, P_pos, isPositive, Y_vj, Dur, ...
   'VariableNames', { ...
     'captainID','Nvoy','productID', ...
     'X1','X2', ...
     'a_c','omega_vj','log_wvj','P_pos','isPositive','Y_vj', 'Duration' } );


end