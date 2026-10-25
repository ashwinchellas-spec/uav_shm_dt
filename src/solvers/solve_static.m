function [d, R] = solve_static(model, F)
%SOLVE_STATIC  Linear static solution  K d = F  with homogeneous Dirichlet BCs.
%   [d, R] = solve_static(model, F)
%   d : full displacement vector (zeros at fixed DOFs)
%   R : reaction forces  R = K d - F  (non-zero only at fixed DOFs)
%   F may have several columns (load cases).
f = model.free;
d = zeros(model.ndof, size(F,2));
d(f,:) = model.K(f,f) \ F(f,:);
if nargout > 1, R = model.K * d - F; end
end
