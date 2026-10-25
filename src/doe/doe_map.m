function [P, names] = doe_map(spec, U)
%DOE_MAP  Map unit-hypercube samples U (N x k) to physical values P (N x k).
%   lin:  x = lo + u (hi - lo)      log:  x = 10^(lo + u (hi - lo))
P = zeros(size(U));
for j = 1:spec.k
    v = spec.vars(j);
    x = v.lo + U(:,j) * (v.hi - v.lo);
    if strcmp(v.scale, 'log'), x = 10.^x; end
    P(:,j) = x;
end
names = {spec.vars.name};
end
