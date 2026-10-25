function [theta, hist] = train_pinet(theta, arch, Xtr, Ttr, Xva, Tva, phys, w, opts)
%TRAIN_PINET  Mini-batch Adam training with early stopping.
%
%   [theta, hist] = train_pinet(theta, arch, Xtr, Ttr, Xva, Tva, phys, w, opts)
%   opts.epochs (300), opts.batch (256), opts.lr (2e-3, cosine-decayed to
%   lr/20), opts.patience (40 epochs), opts.seed, opts.verbose,
%   opts.noise (0): std of Gaussian noise added to the standardised inputs
%   of every mini-batch (data augmentation: the inputs ARE noisy
%   measurements, so the network must not rely on exact values).
%   Model selection: the parameters with the lowest VALIDATION objective
%   are returned.  The validation objective is the supervised loss plus the
%   physics loss with the same weights (the physics loss needs no labels).
%   hist.train / hist.val : objective per epoch;  hist.best_epoch.
def = struct('epochs', 300, 'batch', 256, 'lr', 2e-3, 'patience', 40, 'seed', 1, 'verbose', false, 'noise', 0);
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end, end
rand('state', opts.seed); randn('state', opts.seed); %#ok<RAND>
n = size(Xtr,1);
st = struct('m', 0, 'v', 0, 't', 0);
best = inf; best_theta = theta; hist.train = []; hist.val = []; wait = 0;
fld = fieldnames(Ttr);
for ep = 1:opts.epochs
    lr = opts.lr*(0.05 + 0.95*0.5*(1 + cos(pi*(ep-1)/opts.epochs)));
    perm = randperm(n); Ltot = 0;
    for b = 1:opts.batch:n
        ii = perm(b:min(b+opts.batch-1, n));
        Tb = struct();
        for f = 1:numel(fld), Tb.(fld{f}) = Ttr.(fld{f})(ii,:); end
        Xb = Xtr(ii,:);
        if opts.noise > 0, Xb = Xb + opts.noise*randn(size(Xb)); end
        [Lb, g] = pinet_lossgrad(theta, arch, Xb, Tb, phys, w);
        [theta, st] = adam_step(theta, g, st, lr);
        Ltot = Ltot + Lb*numel(ii);
    end
    Lva = pinet_lossgrad(theta, arch, Xva, Tva, phys, w);
    hist.train(end+1) = Ltot/n; hist.val(end+1) = Lva;
    if Lva < best - 1e-6, best = Lva; best_theta = theta; hist.best_epoch = ep; wait = 0;
    else, wait = wait + 1; end
    if opts.verbose && mod(ep, 25) == 0
        fprintf('    epoch %3d  train %.4f  val %.4f  (best %.4f @ %d)\n', ep, hist.train(end), Lva, best, hist.best_epoch);
    end
    if wait >= opts.patience, break; end
end
theta = best_theta;
end
