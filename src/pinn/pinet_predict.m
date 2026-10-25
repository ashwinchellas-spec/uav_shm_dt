function out = pinet_predict(theta, arch, X, phys, E)
%PINET_PREDICT  Forward pass only: class probabilities, face, damage map,
%   extent, location, reserve factor and decoded features for every row.
%   E : evidence map (n x ne) of evidence_map.m (zeros = no prior).
n = size(X,1);
T = struct('y', ones(n,1), 'lab', false(n,1), 'face', nan(n,1), 'qtrue', zeros(n, arch.ne), ...
           'loc', nan(n,2), 'ext', zeros(n,1), 'rf', zeros(n,1), 'dXn', zeros(n, size(phys.Jn,1)), 'E', E);
w = struct('cls',0,'face',0,'map',0,'loc',0,'ext',0,'rf',0,'phys',0,'l2',0);
[~, ~, ~, out] = pinet_lossgrad(theta, arch, X, T, phys, w);
end
