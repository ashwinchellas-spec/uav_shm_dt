function R = val09_newmark()
%VAL09_NEWMARK  Transient solver checks on a coarse model of the wing.
%   (a) undamped free vibration from the static tip-load shape: the
%       average-acceleration Newmark scheme conserves the total energy
%       exactly for linear systems -> relative energy drift ~ 1e-12.
%   (b) suddenly applied tip load: compare with the exact solution by full
%       modal superposition  d(t) = sum phi_i phi_i' F / w_i^2 (1 - cos w_i t)
%       (dynamic amplification factor -> 2 for an undamped step).
cfg = wing_baseline_config(); cfg.fe.nx = 12; cfg.fe.ny = 3;
model = build_wing_model(cfg);
mesh = model.mesh;
tipn = mesh.sets.tip(round(end/2));
F = point_load(mesh, [cfg.geom.span cfg.geom.c_root/2], 3, 50);
ds = solve_static(model, F);
[fr, Phi, om] = solve_modal(model, numel(model.free));
T1 = 1/fr(1); dt = T1/400; t = 0:dt:3*T1;
wd = 5*(tipn-1) + 3;

o = solve_newmark(model, zeros(model.ndof,1), @(t) zeros(size(t)), t, ...
    struct('d0', ds, 'out_dofs', wd));
drift = max(abs(o.energy/o.energy(1) - 1));
R = vres('V09 Newmark', 'energy drift, undamped free vibration', drift, 0, 1e-9);
R.err = drift; R.pass = drift < R.tol;

o = solve_newmark(model, F, @(t) ones(size(t)), t, struct('out_dofs', wd));
q = (Phi' * F) ./ om.^2;
wex = (Phi(wd,:) .* q') * (1 - cos(om * t));
err = max(abs(o.d - wex)) / max(abs(wex));
R = [R, vres('V09 Newmark', 'step response vs modal superposition (max err)', err, 0, 1e-2)];
R(end).err = err; R(end).pass = err < R(end).tol;
R = [R, vres('V09 Newmark', 'dynamic amplification factor (undamped step)', max(o.d)/ds(wd), 2, 2e-2)];

fig = new_fig(800, 420);
plot(t*1e3, wex*1e3, 'k-', t*1e3, o.d*1e3, 'r--', 'LineWidth', 1.4); grid on;
xlabel('time [ms]'); ylabel('tip deflection [mm]');
legend('exact modal superposition', 'Newmark-\beta (\gamma=1/2, \beta=1/4)', 'Location', 'southeast');
title(sprintf('V09 step response, coarse wing, dt = T_1/400, max error %.2e', err));
save_fig(fig, 'V09_newmark_step');
end
