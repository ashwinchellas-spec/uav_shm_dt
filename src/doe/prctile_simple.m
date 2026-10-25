function p = prctile_simple(x, q)
%PRCTILE_SIMPLE  q-th percentile by linear interpolation (no toolbox needed).
x = sort(x(:)); n = numel(x);
if n == 0, p = nan; return; end
r = 1 + (n - 1)*q/100;
lo = floor(r); hi = ceil(r);
p = x(lo) + (r - lo)*(x(hi) - x(lo));
end
