function [nodes2, weights2] = HG2D(n)
    [x, w] = HG(n);
    [X1, X2] = ndgrid(x, x);
    [W1, W2] = ndgrid(w, w);
    nodes2   = [X1(:), X2(:)];
    weights2 = W1(:) .* W2(:);
end
