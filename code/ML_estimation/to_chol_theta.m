function theta_chol = to_chol_theta(theta) 
% replaces the Sigma matrix (9 elements) by the Cholesky Lower triangular
% matrix, (since theta is positive definite has a Ch. decomposition) 

%  theta: [β(3); α(3); δ(3); γ0; γ1; vec(Sigma_omega) (9)]
Sigma = reshape(theta(12:20), 3, 3); 

if ~isequal(Sigma, Sigma')
    disp('Sigma not symmetric');
  end

L = chol(Sigma, 'lower'); 
theta_chol = [ 
        theta(1:11);      % β1–3, α1–3, δ1–3, γ0, γ1
        L(1,1);       % σ11
        L(2,2);       % σ22
        L(3,3);       % σ33
        L(2,1);       % σ12 = σ21
        L(3,1);       % σ13 = σ31
        L(3,2);         % σ23 = σ32
        theta(21:end) % lambda for each product
    ];
end 