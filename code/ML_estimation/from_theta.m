
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

