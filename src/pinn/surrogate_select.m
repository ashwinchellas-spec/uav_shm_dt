function [Yp, info] = surrogate_select(m, Xr, Xc, Ylog, Plog, tr, va, sopt, seed)
%SURROGATE_SELECT  Train one design-surrogate variant with validation-based L2.
%   [Yp, info] = surrogate_select(m, Xr, Xc, Ylog, Plog, tr, va, sopt, seed)
%   Xr : encoded design variables;  Xc : [Xr, CLT groups of beam_prior]
% m = 2 pure NN, 3 NN + CLT groups, 4 hybrid.  The L2 weight is chosen on
% the VALIDATION designs from {1e-3, 1e-2} (test designs are never used).
best = inf;
for l2 = [1e-3 1e-2]
    o = sopt; o.l2 = l2; o.seed = seed;
    switch m
        case 2, [Y1, i1] = fit_surrogate(Xr, Ylog, Plog, tr, va, 'nn', o);
        case 3, [Y1, i1] = fit_surrogate(Xc, Ylog, Plog, tr, va, 'nn', o);
        case 4, [Y1, i1] = fit_surrogate(Xc, Ylog, Plog, tr, va, 'hybrid', o);
    end
    ev = mean(mean((Y1(va,:) - Ylog(va,:)).^2));
    if ev < best, best = ev; Yp = Y1; info = i1; info.l2 = l2; end
end
end
