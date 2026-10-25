function [X, info] = lhs_maximin(n, k, seed, ncand, nswap)
%LHS_MAXIMIN  Space-filling Latin Hypercube design on [0,1]^k (no toolbox).
%
%   [X, info] = lhs_maximin(n, k, seed, ncand, nswap)
%   n samples, k variables.  Latin Hypercube (McKay et al. 1979): each
%   variable's range is cut into n equal strata and every stratum is used
%   exactly once -> perfect 1-D marginal coverage with any n.
%   Space filling (Morris & Mitchell 1995): of 'ncand' random LHS designs
%   the one with the largest minimum inter-point distance is kept, then
%   improved by 'nswap' greedy column-wise element exchanges that keep the
%   Latin property and are accepted only if the minimum distance grows.
%   Points are jittered uniformly inside their strata.
%
%   info.dmin, info.dmin_random : maximin distance of X and of a plain
%                                 random (Monte-Carlo) design of equal size
if nargin < 3 || isempty(seed), seed = 1; end
if nargin < 4 || isempty(ncand), ncand = 20; end
if nargin < 5 || isempty(nswap), nswap = 200; end
rand('state', seed); %#ok<RAND>   % identical streams in MATLAB and Octave
best = -inf; X = [];
for c = 1:ncand
    P = zeros(n, k);
    for j = 1:k, [~, P(:,j)] = sort(rand(n,1)); end
    Xc = (P - rand(n,k)) / n;
    dm = min_dist(Xc);
    if dm > best, best = dm; X = Xc; end
end
for s = 1:nswap
    j = ceil(rand*k); i = ceil(rand(1,2)*n);
    Y = X; Y(i,j) = Y(fliplr(i),j);
    % jitter stays inside the stratum because both points keep their own
    % stratum index; re-centre within strata
    dm = min_dist(Y);
    if dm > best, best = dm; X = Y; end
end
info.dmin = best;
info.dmin_random = min_dist(rand(n,k));
end

function d = min_dist(X)
n = size(X,1); d = inf;
sq = sum(X.^2, 2);
blk = 500;
for a = 1:blk:n
    ia = a:min(a+blk-1, n);
    D2 = sq(ia) + sq' - 2*X(ia,:)*X';
    for r = 1:numel(ia), D2(r, ia(r)) = inf; end
    d = min(d, sqrt(max(min(D2(:)), 0)));
end
end
