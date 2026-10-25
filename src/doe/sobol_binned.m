function S1 = sobol_binned(X, y, nb)
%SOBOL_BINNED  First-order Sobol indices from a given (LHS) sample.
%   S1 = sobol_binned(X, y, nb)
%   S1_i = Var( E[y | x_i] ) / Var(y)
%   E[y | x_i] is estimated by the mean of y in nb equal-count bins of x_i
%   (no extra model runs; adequate for N >~ 20 nb).  Sum(S1) < 1 means
%   interactions; S1_i ~ 0 means x_i alone does not explain y.
%   Bias correction: with finite n the bin means carry sampling noise, so
%   the raw ratio is biased upward by about b = (nb - 1)/(n - 1) (the
%   expected value for a variable with no influence).  Returned value:
%   S1 = (S_raw - b) / (1 - b), clipped at 0.
if nargin < 3, nb = 10; end
[n, k] = size(X); S1 = zeros(1, k);
vy = var(y, 1);
for i = 1:k
    [~, o] = sort(X(:,i));
    edges = round(linspace(0, n, nb+1));
    m = zeros(nb,1); w = m;
    for b = 1:nb
        ii = o(edges(b)+1:edges(b+1));
        m(b) = mean(y(ii)); w(b) = numel(ii)/n;
    end
    S1(i) = sum(w .* (m - mean(y)).^2) / vy;
end
b = (nb - 1)/(n - 1);
S1 = max((S1 - b)/(1 - b), 0);
end
