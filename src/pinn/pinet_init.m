function [theta, arch] = pinet_init(nin, hidden, ncls, ne, seed)
%PINET_INIT  Parameters of the physics-informed multi-task damage network.
%
%   Trunk (shared):  x -> tanh(W1 x + b1) -> tanh(W2 h1 + b2) = h
%   Heads on h:
%     cls  : ncls logits     -> softmax  p_c   (healthy, 5 damage types, root)
%     face : 1 logit         -> sigmoid  p_top (damaged skin: top vs bottom)
%     map  : ne logits       -> softmax  q_e   (where the damage is; sums to 1)
%     ext  : 1               -> softplus s     (damage extent = sum_e D_e / s_scale)
%     rf   : 1               -> linear         (residual strength, standardised)
%   Map prior: the map logits receive  + a*E  where E is a physics-based
%   evidence map (evidence_map.m) and a is a learnable gain (initial 1).
%   Glorot-uniform initialisation; parameter cells are packed into one
%   vector theta (nn_pack) so Adam and the gradient check are generic.
rand('state', seed); %#ok<RAND>
dims = [nin, hidden(:)'];
P = {};
for l = 1:numel(hidden)
    a = sqrt(6/(dims(l) + dims(l+1)));
    P{end+1} = a*(2*rand(dims(l), dims(l+1)) - 1); %#ok<AGROW>
    P{end+1} = zeros(1, dims(l+1)); %#ok<AGROW>
end
h = dims(end);
heads = [ncls, 1, ne, 1, 1];
for k = 1:numel(heads)
    a = sqrt(6/(h + heads(k)));
    P{end+1} = a*(2*rand(h, heads(k)) - 1); %#ok<AGROW>
    P{end+1} = zeros(1, heads(k)); %#ok<AGROW>
end
P{end+1} = 1;                                       % gain a of the evidence-map prior
[theta, arch.shapes] = nn_pack(P);
arch.nin = nin; arch.hidden = hidden; arch.ncls = ncls; arch.ne = ne;
arch.nlayer = numel(hidden);
end
