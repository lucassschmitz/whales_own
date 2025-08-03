function [x, w] = IE9_hermiteGaussRule(n)
    i = (1:n-1)';
    a = sqrt(i/2);
    T = diag(a,1) + diag(a,-1);
    [V, D] = eig(T);
    [x_sorted, idx] = sort(diag(D));     % GH nodes
    V = V(:, idx);
    x = x_sorted;
    w = (sqrt(pi) * (V(1,:).^2))';       % ensure w is column (n×1)
end