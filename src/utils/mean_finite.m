function m = mean_finite(v)
%MEAN_FINITE  Mean of the finite entries of v (toolbox-free nanmean).
v = v(isfinite(v));
if isempty(v), m = NaN; else, m = mean(v); end
end
