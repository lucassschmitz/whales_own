function [nodes, weights] = HG3D(n)
    %returns an (n^3×3) array of 3D nodes and an (n^3×1) vector of weights

    [x, w] = HG(n);

    % Build 3‐D grids of nodes and weights
    [X1, X2, X3] = ndgrid(x, x, x);
    [W1, W2, W3] = ndgrid(w, w, w);

    % Pack into list of 3D nodes and corresponding weights
    nodes   = [X1(:), X2(:), X3(:)];
    weights = W1(:) .* W2(:) .* W3(:);
end
