
function Tsim = simulate_partial_data(T, J, s_omega, alpha, delta, beta, gamma0, gamma1, lambda)

%differences with IE11_gen_data: receives input lambda and uses to multiply
%production. 
%1.instead of simulating the durations and the tonnage uses the actual
%duration and tonnage. Receives as input T which is the table with the
%actual data. 

 % beta      - J×2 matrix of covariate coefficients (product-specific)  <-- CHANGED


%find unique captains
capt_vec = unique(T.captainID);
C = numel(capt_vec);
N  = height(T);
voyages = unique(T.voyageID); 

%Pre‐allocate the “long” vectors exactly as before
captainID = zeros(N,1);
Nvoy = zeros(N,1);
productID = zeros(N,1);

X1        = zeros(N,1);

a_c       = zeros(N,1);
omega_vj  = zeros(N,1);   % STILL a scalar per row, but filled from a J‐vector
log_wvj   = zeros(N,1);
P_pos     = zeros(N,1);
isPositive= zeros(N,1);
Y_vj      = zeros(N,1);
Dur       = zeros(N,1); 

idx = 0;


%Pre‐draw all a_c for each captain
true_ac = randn(C,1);


for k = 1:numel(voyages)
    vid = voyages(k); % voyage id 
    rows = find(T.voyageID == vid);
    %Draw J‐vector of omegas for this voyage
    Omega_vec = mvnrnd(zeros(1,J), s_omega');
    
    %Fill in each product‐row for this voyage
    for r = rows'
        c   = T.captainID(r);
        j   = T.productID(r);
        x1  = T.X1(r);
        tau = T.Tau(r);
        
        captainID(r) = c;
        Nvoy(r)      = T.n_voy(r);
        productID(r) = j;
        X1(r)        = x1;
        X2(r)        = rand > 0.5;
        a_c(r)       = true_ac(c);
        omega_vj(r)  = Omega_vec(j);
        
        %Linear index utility
        log_wvj(r) = beta(j,1)*x1 + delta(j)*a_c(r) + omega_vj(r);
        %Probability of positive outcome
        P_pos(r) = 1/(1 + exp(gamma0 - gamma1 * log_wvj(r)));
        isPositive(r) = rand < P_pos(r);
        
        %Generate output Y
        if isPositive(r)
            Y_vj(r) = lambda(j)* (tau * exp(log_wvj(r)))^alpha(j);
        else
            Y_vj(r) = 0;
        end
        Dur(r) = tau;
    end
end

%Pack results into a table
Tsim = table( captainID, Nvoy, productID, X1, a_c, omega_vj, ...
               log_wvj, P_pos, isPositive, Y_vj, Dur, ...
               'VariableNames', { 'captainID','Nvoy','productID', ...
                                   'X1','a_c','omega_vj', ...
                                   'log_wvj','P_pos','isPositive', ...
                                   'Y_vj','Duration' } );

end
