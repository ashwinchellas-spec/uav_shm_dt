function [freq, Phi, omega] = solve_modal(model, nm)
%SOLVE_MODAL  Natural frequencies and mass-normalised mode shapes.
%
%   [freq, Phi, omega] = solve_modal(model, nm)
%   Solves  (K - omega^2 M) phi = 0  on the free DOFs.
%   freq  [Hz] ascending, omega [rad/s],
%   Phi   ndof x nm, phi' M phi = 1, sign fixed so that the largest |w|
%         component is positive (removes the arbitrary eigenvector sign,
%         which matters for mode-shape features / ML in later milestones).
%
%   Robustness: the iterative solver (eigs, shift-invert about 0) can fail
%   to converge for some eigenvalues; MATLAB then returns NaN WITHOUT an
%   error.  Every result is therefore checked (finite, positive, nm values);
%   on failure eigs is retried with a larger Krylov subspace and a looser
%   tolerance, and finally the exact dense solver eig() is used.

f = model.free;
Kf = model.K(f,f); Mf = model.M(f,f);
n = numel(f);
ok = false;
if n > 1500
    attempts = {struct('tol', 1e-12, 'maxit', 1000, 'p', min(n, max(2*nm + 1, 20))), ...
                struct('tol', 1e-10, 'maxit', 3000, 'p', min(n, max(4*nm, 60)))};
    for a = 1:numel(attempts)
        opts = attempts{a}; opts.disp = 0;
        try
            ws = warning('off', 'all');
            [V, L] = eigs(Kf, Mf, nm, 0, opts);
            warning(ws);
            lam = real(diag(L));
            ok = numel(lam) == nm && all(isfinite(lam)) && all(lam > 0) && all(isfinite(V(:)));
        catch
            warning(ws); ok = false;
        end
        if ok, break; end
    end
end
if ~ok
    [V, L] = eig(full(Kf), full(Mf));
    lam = diag(L);
end
lam = real(lam); V = real(V);
[lam, i] = sort(lam); V = V(:, i);
lam = lam(1:nm); V = V(:, 1:nm);

Phi = zeros(model.ndof, nm);
Phi(f,:) = V;
wdof = 3:5:model.ndof;
for k = 1:nm
    Phi(:,k) = Phi(:,k) / sqrt(Phi(:,k)' * model.M * Phi(:,k));
    [~, j] = max(abs(Phi(wdof,k)));
    if Phi(wdof(j),k) < 0, Phi(:,k) = -Phi(:,k); end
end
omega = sqrt(max(lam, 0));
freq  = omega / (2*pi);
end
