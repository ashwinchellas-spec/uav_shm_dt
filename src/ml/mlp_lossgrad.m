function [L, g, Yp] = mlp_lossgrad(theta, arch, X, Y, l2)
%MLP_LOSSGRAD  Mean-squared error (+ L2 on weights) and its exact gradient.
%   [L, g, Yp] = mlp_lossgrad(theta, arch, X, Y, l2)
%   Hidden layers tanh, output linear.  Y may be [] (prediction only).
P = nn_unpack(theta, arch.shapes); nl = arch.nlayer;
A = cell(nl+1, 1); A{1} = X;
for l = 1:nl
    Z = A{l}*P{2*l-1} + P{2*l};
    if l < nl, A{l+1} = tanh(Z); else, A{l+1} = Z; end
end
Yp = A{end};
if isempty(Y), L = []; g = []; return; end
n = size(X,1); R = Yp - Y;
W2 = 0; for l = 1:nl, W2 = W2 + sum(P{2*l-1}(:).^2); end
L = sum(R(:).^2)/(n*size(Y,2)) + l2*W2;
if nargout < 2, return; end
dP = cell(size(P));
dZ = 2*R/(n*size(Y,2));
for l = nl:-1:1
    dP{2*l-1} = A{l}'*dZ + 2*l2*P{2*l-1};
    dP{2*l} = sum(dZ, 1);
    if l > 1, dZ = (dZ*P{2*l-1}') .* (1 - A{l}.^2); end
end
g = nn_pack(dP);
end
