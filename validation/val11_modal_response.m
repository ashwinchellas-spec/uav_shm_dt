function R = val11_modal_response()
%VAL11_MODAL_RESPONSE  Fast modal solver (used for every DOE sample) vs Newmark.
%   Damped 1-cos gust response of the baseline wing (t_g = T1 to excite
%   dynamics), tip deflection and acceleration:
%     reference : full-model Newmark-beta, dt = T1/2000 (converged)
%     test      : 8-mode exact-FOH modal superposition + mode-acceleration
%                 static correction, dt = T1/100
cfg = wing_baseline_config(); model = build_wing_model(cfg);
[fr, Phi, om] = solve_modal(model, 8);
[al, be] = rayleigh_coeffs(om(1), om(5), 0.02);
zeta = al./(2*om) + be*om/2;
[lc] = wing_load_cases(cfg, model); F0 = lc(4).F;
tip = model.mesh.sets.tip; [~, j] = min(abs(model.mesh.eta(tip) - 0.5));
wd = 5*(tip(j)-1) + 3;
T1 = 1/fr(1); tg = T1;
gf = @(t) 0.5*(1 - cos(2*pi*t/tg)) .* (t <= tg);

t1 = 0:T1/2000:4*T1;
o1 = solve_newmark(model, F0, gf, t1, struct('alpha', al, 'beta_r', be, 'out_dofs', wd));
t2 = 0:T1/100:4*T1;
ds = solve_static(model, F0);
gam = Phi' * F0;
resid = ds(wd) - Phi(wd,:) * (gam ./ om.^2);
tic; o2 = modal_response(om, zeta, Phi(wd,:), gam, gf(t2), t2(2)-t2(1), resid); tm = toc;
w1 = interp1(t1, o1.d, t2); a1 = interp1(t1, o1.a, t2);
ed = max(abs(o2.d - w1)) / max(abs(w1));
ea = max(abs(o2.a - a1)) / max(abs(a1));
R = [vres('V11 modal response', 'tip deflection, 8-mode FOH vs Newmark (max err)', ed, 0, 1e-2), ...
     vres('V11 modal response', 'tip acceleration, 8-mode FOH vs Newmark (max err)', ea, 0, 5e-2), ...
     vres('V11 modal response', sprintf('modal solver time [s] (%d steps)', numel(t2)), tm, 0, inf)];
R(1).err = ed; R(1).pass = ed < R(1).tol; R(2).err = ea; R(2).pass = ea < R(2).tol;
R(3).err = tm;
fig = new_fig(900, 600);
subplot(2,1,1); plot(t1*1e3, o1.d*1e3, 'k-', t2*1e3, o2.d*1e3, 'r.', 'LineWidth', 1.2); grid on;
ylabel('tip w [mm]'); legend('Newmark full model', '8-mode FOH + static correction', 'Location', 'northeast');
title(sprintf('V11: gust t_g = T_1, max errors: w %.1e, a %.1e', ed, ea));
subplot(2,1,2); plot(t1*1e3, o1.a, 'k-', t2*1e3, o2.a, 'r.', 'LineWidth', 1.2); grid on;
xlabel('time [ms]'); ylabel('tip acceleration [m/s^2]');
save_fig(fig, 'V11_modal_response');
end
