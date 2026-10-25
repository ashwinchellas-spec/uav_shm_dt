function [L, g, parts, out] = pinet_lossgrad(theta, arch, X, T, phys, w)
%PINET_LOSSGRAD  Loss and exact gradient of the physics-informed damage network.
%
%   [L, g, parts, out] = pinet_lossgrad(theta, arch, X, T, phys, w)
%   X (n x nin) standardised features.  T: targets (see below).  phys:
%   linear FE decoder.  w: loss weights.  Backpropagation is derived by
%   hand (no toolbox); it is verified by grad_check to ~1e-8.
%
%   SUPERVISED terms (only rows with T.lab = true):
%     L_cls  = -mean log p_c(y)                              class
%     L_face = BCE(p_top, face)            (skin / delam classes)
%     L_map  = -mean sum_e q*_e log q_e    (damage-map shape, q* = D/sum D)
%     L_loc  = mean |sum_e q_e x_e - x*|^2 / ell^2           centre
%     L_ext  = mean (s - s*)^2             (extent  s* = sum_e D_e / s_scale)
%     L_rf   = mean (rf - rf*)^2           (standardised reserve factor)
%   PHYSICS term (ALL rows, needs no labels):
%     predicted damage field   d_e = s_scale * s * q_e
%     decoded feature change   x^ = sum_t pi_t J_t d  +  pi_root J_root s_scale s
%     L_phys = mean_n mean_k ((x^ - dx)/sigma)^2
%     pi_t = type probabilities from p_c and p_top (skin damage on the top
%     or bottom face changes different sensors).  J_t: FE sensitivities
%     (damage_sensitivities.m) -> the network's diagnosis must be able to
%     REPRODUCE what the sensors measured, through the governing equations.
%
%   T fields: y, lab, face (NaN = n/a), qtrue (n x ne), loc (n x 2, NaN =
%   n/a), ext, rf, dXn (n x nf, measured feature change / sigma), E (n x ne
%   physics evidence map: map logits = W h + b + a E, a learnable).
P = nn_unpack(theta, arch.shapes);
n = size(X,1); nl = arch.nlayer;
% ---- forward --------------------------------------------------------------
A = cell(nl+1,1); A{1} = X;
for l = 1:nl
    A{l+1} = tanh(A{l}*P{2*l-1} + P{2*l});
end
H = A{end}; k0 = 2*nl;
Lc = H*P{k0+1} + P{k0+2};  Lf = H*P{k0+3} + P{k0+4};  Lm = H*P{k0+5} + P{k0+6} + P{end}*T.E;
Le = H*P{k0+7} + P{k0+8};  Lr = H*P{k0+9} + P{k0+10};
pc = exp(Lc - max(Lc,[],2)); pc = pc ./ sum(pc,2);
pf = 1 ./ (1 + exp(-Lf));
q  = exp(Lm - max(Lm,[],2)); q = q ./ sum(q,2);
s  = log(1 + exp(-abs(Le))) + max(Le, 0);           % stable softplus
sg = 1 ./ (1 + exp(-Le));                            % d softplus / d Le
rf = Lr;
locp = q * phys.XY;

% ---- losses --------------------------------------------------------------------
lab = logical(T.lab(:)); nlb = max(sum(lab), 1);
Y = zeros(n, arch.ncls); Y(sub2ind(size(Y), (1:n)', T.y(:))) = 1;
parts.cls = -sum(log(sum(pc(lab,:).*Y(lab,:), 2) + 1e-12)) / nlb;
fm = lab & isfinite(T.face); nfm = max(sum(fm), 1);
parts.face = -sum(T.face(fm).*log(pf(fm) + 1e-12) + (1 - T.face(fm)).*log(1 - pf(fm) + 1e-12)) / nfm;
mm = lab & T.y >= 2 & T.y <= 6; nmm = max(sum(mm), 1);
parts.map = -sum(sum(T.qtrue(mm,:).*log(q(mm,:) + 1e-12))) / nmm;
lm = lab & all(isfinite(T.loc), 2); nlm = max(sum(lm), 1);
dl = locp(lm,:) - T.loc(lm,:);
parts.loc = sum(dl(:).^2) / (phys.ell^2 * nlm);
parts.ext = sum((s(lab) - T.ext(lab)).^2) / nlb;
parts.rf  = sum((rf(lab) - T.rf(lab)).^2) / nlb;
% physics decoder
pi_ = [pc(:,2), pc(:,3).*(1-pf), pc(:,3).*pf, pc(:,4).*(1-pf), pc(:,4).*pf, pc(:,5), ...
       pc(:,6).*(1-pf), pc(:,6).*pf, pc(:,7)];
dh = phys.s_scale * s .* q;                          % predicted damage field (n x ne)
nt = size(phys.Jn, 3);
C = cell(nt+1, 1); xh = zeros(n, size(phys.Jn,1));
for t = 1:nt
    C{t} = dh * phys.Jn(:,:,t)';
    xh = xh + pi_(:,t) .* C{t};
end
C{nt+1} = (phys.s_scale * s) * phys.Jrootn';
xh = xh + pi_(:,nt+1) .* C{nt+1};
r = xh - T.dXn;
nfp = size(r, 2);
parts.phys = sum(r(:).^2) / (n * nfp);
W2sum = 0;
for l = 1:nl, W2sum = W2sum + sum(P{2*l-1}(:).^2); end
for k = 1:5,  W2sum = W2sum + sum(P{k0+2*k-1}(:).^2); end
parts.l2 = W2sum;
L = w.cls*parts.cls + w.face*parts.face + w.map*parts.map + w.loc*parts.loc + ...
    w.ext*parts.ext + w.rf*parts.rf + w.phys*parts.phys + w.l2*parts.l2;

if nargout > 3
    out = struct('pc', pc, 'pf', pf, 'q', q, 's', s, 'rf', rf, 'loc', locp, 'xh', xh);
end
if nargout < 2, return; end

% ---- backward ------------------------------------------------------------------------
G = 2*w.phys*r / (n*nfp);                            % dL/dxh
dpi = zeros(n, nt+1); ddh = zeros(n, arch.ne); ds = zeros(n,1);
for t = 1:nt
    dpi(:,t) = sum(G .* C{t}, 2);
    ddh = ddh + pi_(:,t) .* (G * phys.Jn(:,:,t));
end
dpi(:,nt+1) = sum(G .* C{nt+1}, 2);
ds = ds + phys.s_scale * pi_(:,nt+1) .* (G * phys.Jrootn);
dq = phys.s_scale * s .* ddh;
ds = ds + phys.s_scale * sum(ddh .* q, 2);
% type probabilities -> p_c, p_top
dpc = zeros(n, arch.ncls); dpf = zeros(n,1);
dpc(:,2) = dpi(:,1);
dpc(:,3) = dpi(:,2).*(1-pf) + dpi(:,3).*pf;   dpf = dpf + pc(:,3).*(dpi(:,3) - dpi(:,2));
dpc(:,4) = dpi(:,4).*(1-pf) + dpi(:,5).*pf;   dpf = dpf + pc(:,4).*(dpi(:,5) - dpi(:,4));
dpc(:,5) = dpi(:,6);
dpc(:,6) = dpi(:,7).*(1-pf) + dpi(:,8).*pf;   dpf = dpf + pc(:,6).*(dpi(:,8) - dpi(:,7));
dpc(:,7) = dpi(:,9);
dLc = pc .* (dpc - sum(dpc .* pc, 2));
dLc(lab,:) = dLc(lab,:) + w.cls*(pc(lab,:) - Y(lab,:)) / nlb;
dLf = dpf .* pf .* (1 - pf);
dLf(fm) = dLf(fm) + w.face*(pf(fm) - T.face(fm)) / nfm;
% map: CE on shape + location of the centroid
dq(lm,:) = dq(lm,:) + w.loc * 2*dl*phys.XY' / (phys.ell^2 * nlm);
dLm = q .* (dq - sum(dq .* q, 2));
dP = cell(size(P));
dLm(mm,:) = dLm(mm,:) + w.map*(q(mm,:) - T.qtrue(mm,:)) / nmm;
dP{numel(P)} = sum(sum(dLm .* T.E));
ds(lab) = ds(lab) + w.ext*2*(s(lab) - T.ext(lab)) / nlb;
dLe = ds .* sg;
dLr = zeros(n,1); dLr(lab) = w.rf*2*(rf(lab) - T.rf(lab)) / nlb;
% heads
hd = {dLc, dLf, dLm, dLe, dLr};
dH = zeros(size(H));
for k = 1:5
    dP{k0+2*k-1} = H' * hd{k} + 2*w.l2*P{k0+2*k-1};
    dP{k0+2*k}   = sum(hd{k}, 1);
    dH = dH + hd{k} * P{k0+2*k-1}';
end
% trunk
dA = dH;
for l = nl:-1:1
    dZ = dA .* (1 - A{l+1}.^2);
    dP{2*l-1} = A{l}' * dZ + 2*w.l2*P{2*l-1};
    dP{2*l}   = sum(dZ, 1);
    dA = dZ * P{2*l-1}';
end
g = nn_pack(dP);
end
