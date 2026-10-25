function [theta, hist] = train_mlp(theta, arch, Xtr, Ytr, Xva, Yva, opts)
%TRAIN_MLP  Full-batch / mini-batch Adam for mlp_lossgrad with early stopping.
%   opts.epochs (3000), opts.batch (64), opts.lr (3e-3, cosine to lr/20),
%   opts.l2 (1e-4), opts.patience (300), opts.seed (1).
%   Returns the parameters with the lowest validation MSE.
def = struct('epochs', 3000, 'batch', 64, 'lr', 3e-3, 'l2', 1e-4, 'patience', 300, 'seed', 1);
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end, end
rand('state', opts.seed); %#ok<RAND>
n = size(Xtr,1); st = struct('m', 0, 'v', 0, 't', 0);
best = inf; best_theta = theta; wait = 0; hist.val = []; hist.best_epoch = 0;
for ep = 1:opts.epochs
    lr = opts.lr*(0.05 + 0.95*0.5*(1 + cos(pi*(ep-1)/opts.epochs)));
    perm = randperm(n);
    for b = 1:opts.batch:n
        ii = perm(b:min(b+opts.batch-1, n));
        [~, g] = mlp_lossgrad(theta, arch, Xtr(ii,:), Ytr(ii,:), opts.l2);
        [theta, st] = adam_step(theta, g, st, lr);
    end
    [~, ~, Yp] = mlp_lossgrad(theta, arch, Xva, [], 0);
    Lva = mean((Yp(:) - Yva(:)).^2); hist.val(end+1) = Lva;
    if Lva < best - 1e-9, best = Lva; best_theta = theta; hist.best_epoch = ep; wait = 0;
    else, wait = wait + 1; end
    if wait >= opts.patience, break; end
end
theta = best_theta;
end
