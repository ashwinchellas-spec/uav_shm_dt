function [theta, arch] = mlp_init(nin, hidden, nout, seed)
%MLP_INIT  Fully connected regression network  x -> tanh layers -> linear.
%   [theta, arch] = mlp_init(nin, hidden, nout, seed)   (Glorot uniform)
rand('state', seed); %#ok<RAND>
dims = [nin, hidden(:)', nout];
P = {};
for l = 1:numel(dims)-1
    a = sqrt(6/(dims(l) + dims(l+1)));
    P{end+1} = a*(2*rand(dims(l), dims(l+1)) - 1); %#ok<AGROW>
    P{end+1} = zeros(1, dims(l+1)); %#ok<AGROW>
end
[theta, arch.shapes] = nn_pack(P);
arch.nlayer = numel(dims) - 1;
end
