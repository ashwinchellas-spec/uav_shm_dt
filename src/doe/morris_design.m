function [U, delta] = morris_design(k, r, p, seed)
%MORRIS_DESIGN  Morris (1991) elementary-effects trajectories on [0,1]^k.
%   [U, delta] = morris_design(k, r, p, seed)
%   r trajectories of k+1 points each (r(k+1) model runs), p grid levels,
%   step delta = p / (2(p-1)).  In every trajectory each factor moves
%   exactly once by +/- delta, in random order, from a random grid start.
if nargin < 3, p = 4; end
if nargin < 4, seed = 7; end
rand('state', seed); %#ok<RAND>
delta = p/(2*(p-1));
U = zeros(r*(k+1), k);
B = tril(ones(k+1, k), -1);
for t = 1:r
    xs = floor(rand(1,k)*(p/2)) / (p-1);        % start levels allowing +delta
    D = diag(2*(rand(1,k) > 0.5) - 1);
    [~, perm] = sort(rand(1,k)); Pm = eye(k); Pm = Pm(:, perm);
    J = ones(k+1, 1);
    Bs = (J*xs + (delta/2)*((2*B - J*ones(1,k))*D + J*ones(1,k))) * Pm;
    U((t-1)*(k+1) + (1:k+1), :) = Bs;
end
end
