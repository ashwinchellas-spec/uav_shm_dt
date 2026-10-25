function rel = grad_check(fun, theta, nchk, h)
%GRAD_CHECK  Compare analytic gradient with central finite differences.
%   rel = grad_check(fun, theta, nchk, h)
%   fun : handle [L, g] = fun(theta)
%   Checks nchk randomly chosen components; returns the maximum relative
%   error |g_a - g_fd| / max(|g_a| + |g_fd|, 1e-8).  A correct backprop
%   gives ~1e-7 or smaller (double precision, h = 1e-6).
if nargin < 3, nchk = 30; end
if nargin < 4, h = 1e-6; end
[~, g] = fun(theta);
idx = randperm(numel(theta), min(nchk, numel(theta)));
rel = 0;
for i = idx
    tp = theta; tp(i) = tp(i) + h; Lp = fun(tp);
    tm = theta; tm(i) = tm(i) - h; Lm = fun(tm);
    gfd = (Lp - Lm) / (2*h);
    rel = max(rel, abs(g(i) - gfd) / max(abs(g(i)) + abs(gfd), 1e-8));
end
end
