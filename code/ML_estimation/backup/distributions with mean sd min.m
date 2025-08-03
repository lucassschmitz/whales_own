clear 
clc

% Number of draws
N = 1000000;

% parameters obtained from our data 

%1. duration (taken as continuous; mean 2.5, s.d 1.5, min 1) 
mu_tau = 2.5; 
sd_tau = 1.5; 
min_tau = 1; 

%2. tonnage (mean 272, sd 103, min 45 )
mu_x = 272;
sd_x   = 103;
min_x  = 15; % set 15 to get around 45 


mu_Z = mu_tau - min_tau; 
sigma2 = log(1 + (sd_tau/mu_Z)^2); 
mu = log(mu_Z) - sigma2/2; 
Tau = lognrnd(mu, sqrt(sigma2), 1, 1); 
Tau = Tau + min_tau;

[min(Tau), mean(Tau), std(Tau) ]

mu_Z = mu_x - min_x; 
sigma2 = log(1 + (sd_x/mu_Z)^2); 
mu = log(mu_Z) - sigma2/2; 
X = lognrnd(mu, sqrt(sigma2), N, 1); 
X = X + min_x; 

[min(X), mean(X), std(X) ]