function theta = to_theta(theta_chol) 

L = [theta_chol(12), 0, 0;  
    theta_chol(15), theta_chol(13), 0; 
    theta_chol(16), theta_chol(17), theta_chol(14)];
sigma = L * L';

theta = [theta_chol(1:11);  % copy β, α, δ, γ0, γ1
        reshape(sigma,9,1); 
        theta_chol(end-2:end)];

end


 % v1 commented below 
% function theta = to_theta(theta_red) 
%     %function to create the full vector based on the reduced vector. 
% 
%     theta =  [    theta_red(1:11);           % copy β, α, δ, γ0, γ1
%     reshape([                   % reconstruct Σ_omega columnwise
%         theta_red(12), theta_red(15), theta_red(16);  ...
%         theta_red(15), theta_red(13), theta_red(17);  ...
%         theta_red(16), theta_red(17), theta_red(14)       
%     ], 9, 1)
%     ];
% 
% end