function [Yp, info] = fit_surrogate(Xd, Ylog, Plog, tr, va, mode, opts)
%FIT_SURROGATE  Pure data-driven or hybrid physics-AI surrogate.
%
%   [Yp, info] = fit_surrogate(Xd, Ylog, Plog, tr, va, mode, opts)
%   Xd   : encoded design variables (design_inputs.m), N x p
%   Ylog : log of the FE responses (N x m);  Plog : log of the reduced-order
%          physics prediction (beam_prior.m), N x m
%   tr, va : training / validation row masks;  opts: train_mlp options +
%          opts.hidden (default [32 32]), opts.seed
%   mode : 'nn'     -  y = NN(x)                       (black box)
%          'hybrid' -  y = y_phys(x) * exp(NN(x, log y_phys))
%                      the network learns only the MODEL-FORM ERROR of the
%                      physics model (residual / discrepancy learning, Kennedy
%                      & O'Hagan 2001); with no data it falls back to physics.
%   Yp   : predicted log responses for all rows
if ~isfield(opts, 'hidden'), opts.hidden = [32 32]; end
if ~isfield(opts, 'seed'), opts.seed = 1; end
switch mode
    case 'nn',     Xin = Xd;          T = Ylog;
    case 'hybrid', Xin = [Xd, Plog];  T = Ylog - Plog;
    otherwise, error('fit_surrogate: unknown mode');
end
mx = mean(Xin(tr,:), 1); sx = std(Xin(tr,:), 0, 1); sx(sx == 0) = 1;
mt = mean(T(tr,:), 1);   st = std(T(tr,:), 0, 1);  st(st == 0) = 1;
Z = (Xin - mx) ./ sx; Tz = (T - mt) ./ st;
[th, arch] = mlp_init(size(Z,2), opts.hidden, size(T,2), opts.seed);
[th, h] = train_mlp(th, arch, Z(tr,:), Tz(tr,:), Z(va,:), Tz(va,:), opts);
[~, ~, Tp] = mlp_lossgrad(th, arch, Z, [], 0);
Yp = Tp.*st + mt;
if strcmp(mode, 'hybrid'), Yp = Yp + Plog; end
info = struct('theta', th, 'arch', arch, 'best_epoch', h.best_epoch, 'mx', mx, 'sx', sx, 'mt', mt, 'st', st);
end
