%RUN_MILESTONE2  Step 2 - in-house FSDT solver applied to the baseline wing.
%
%   (Verification & validation of the solver itself: run_all_validation.)
%   This script produces every wing result quoted in Milestone2_Report.md:
%     A. mesh-convergence study + Richardson extrapolation
%     B. static analysis of all design load cases (deflection, twist,
%        strains, Tsai-Wu / Hashin / core-shear / interlaminar indices,
%        reserve factors, global equilibrium check)
%     C. importance of transverse shear (FSDT vs shear-rigid vs k = 5/6)
%     D. modal analysis (frequencies, mode identification, orthogonality)
%     E. dynamic 1-cos gust response (Newmark) and gust-duration sweep
%     F. damage sensitivity demonstration (bridge to Milestone 3)
%     G. computational cost (for DOE planning)
%   Outputs: results/figures/M2_*.png, results/data/M2_summary.txt,
%            results/data/M2_results.mat
%
%   >> run_milestone2

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
cfg = wing_baseline_config();
fid = fopen(fullfile(root, 'results', 'data', 'M2_summary.txt'), 'w');
P = @(varargin) logprintf(fid, varargin{:});
P('=====================================================================\n');
P(' MILESTONE 2 - FSDT solver results for %s\n', cfg.name);
P('=====================================================================\n');
c = cfg.geom.c_root;

%% A. Mesh convergence ----------------------------------------------------
P('\n[A] MESH CONVERGENCE (LC4 critical gust, tip mid-chord deflection; f1..f6)\n');
ks = [2 4 8 16];                       % nx = 4k, ny = k (aspect-ratio-1 elements)
meshes = [4*ks(:) ks(:); cfg.fe.nx cfg.fe.ny];
conv = zeros(size(meshes,1), 9);
x_off = 0.05;       % stress evaluation station: >= 2h from the clamped root
for i = 1:size(meshes,1)
    cf = cfg; cf.fe.nx = meshes(i,1); cf.fe.ny = meshes(i,2);
    mdl = build_wing_model(cf);
    lci = wing_load_cases(cf, mdl);
    d = solve_static(mdl, lci(4).F);
    tip = mdl.mesh.sets.tip;
    [~, j] = min(abs(mdl.mesh.eta(tip) - 0.5));
    fr = solve_modal(mdl, 6);
    rs = recover_stresses(mdl, d);
    xc = elem_centroids(mdl.mesh);
    tw_e = max(rs.tsaiwu, [], 1);
    % off-corner metric: Tsai-Wu at the station x = x_off, interpolated
    % linearly between element-centroid columns (mesh independent location)
    ux = unique(round(xc(:,1)*1e9)/1e9);
    twcol = arrayfun(@(u) max(tw_e(abs(xc(:,1) - u) < 1e-9)), ux);
    tw_off = interp1(ux, twcol, x_off);
    conv(i,:) = [d(5*(tip(j)-1)+3), fr(:)', max(tw_e), tw_off];
    P('  %3dx%-3d (%5d DOF): w_tip = %.5f mm  f = [%s] Hz  TW_max = %.4f  TW(x=%.0fmm) = %.4f\n', ...
      meshes(i,1), meshes(i,2), mdl.ndof, conv(i,1)*1e3, sprintf('%.3f ', fr), conv(i,8), x_off*1e3, conv(i,9));
end
% Richardson extrapolation on the three finest ratio-2 meshes
q3 = conv(2:4, :);
pobs = log(abs((q3(2,:) - q3(1,:)) ./ (q3(3,:) - q3(2,:)))) / log(2);
pobs(~isfinite(pobs)) = 2;
qinf = q3(3,:) + (q3(3,:) - q3(2,:)) ./ (2.^pobs - 1);
ebase = abs(conv(end,:) ./ qinf - 1) * 100;
names = {'w_tip', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'TW_max', 'TW(x=50mm)'};
P('  Richardson (16x4, 32x8, 64x16):\n');
units = {'mm', 'Hz', 'Hz', 'Hz', 'Hz', 'Hz', 'Hz', '-', '-'};
scl = [1e3 1 1 1 1 1 1 1 1];
for k = 1:9
    P('    %-11s observed order p = %5.2f, extrapolated = %10.5g %-3s baseline 40x10 error = %.3f %%\n', ...
      names{k}, pobs(k), qinf(k)*scl(k), units{k}, ebase(k));
end
P('  NOTE: TW_max sits at the clamped trailing-edge root corner, where a clamped\n');
P('  and a free edge meet: plate theory predicts a stress singularity there, so the\n');
P('  peak grows with refinement (negative observed order). Design values are\n');
P('  therefore read at x = %.0f mm (about 2h from the root), which converges.\n', x_off*1e3);
fig = new_fig(800, 500);
hsz = cfg.geom.span ./ meshes(1:4,1);
loglog(hsz*1e3, abs(conv(1:4,1)/qinf(1) - 1), 'bo-', hsz*1e3, abs(conv(1:4,2)/qinf(2) - 1), 'rs-', ...
       hsz*1e3, abs(conv(1:4,5)/qinf(5) - 1), 'g^-', hsz*1e3, abs(conv(1:4,9)/qinf(9) - 1), 'mv-', 'LineWidth', 1.5);
grid on; xlabel('element size [mm]'); ylabel('relative error vs Richardson limit');
legend('tip deflection LC4', 'f_1', 'f_4', 'Tsai-Wu at x = 50 mm', 'Location', 'southeast');
title('M2-A mesh convergence of the wing model');
save_fig(fig, 'M2_mesh_convergence');

%% Baseline model --------------------------------------------------------
tA = tic; model = build_wing_model(cfg); t_asm = toc(tA);
[lc, env] = wing_load_cases(cfg, model);
mesh = model.mesh; nn = mesh.nn; X = mesh.nodes(:,1); Y = mesh.nodes(:,2);
tip = mesh.sets.tip; [~, order] = sort(mesh.eta(tip)); tip = tip(order);
tipLE = tip(1); tipTE = tip(end); [~, j] = min(abs(mesh.eta(tip) - 0.5)); tipMC = tip(j);
wd = @(n) 5*(n-1) + 3;

%% B. Static analysis ------------------------------------------------------
P('\n[B] STATIC ANALYSIS, mesh %dx%d (%d DOF)\n', mesh.nx, mesh.ny, model.ndof);
P('  %-19s %6s %9s %9s %9s %8s %9s %9s %8s %8s %8s %8s\n', 'case', 'n', 'wLE[mm]', 'wMC[mm]', ...
  'wTE[mm]', 'twist', 'ex_max', 'ex_min', 'TsaiWu', 'Hashin', 'Core-S', 'ILSS');
P('  %-19s %6s %9s %9s %9s %8s %9s %9s %8s %8s %8s %8s\n', '', '', '', '', '', '[deg]', '[ue]', '[ue]', 'index', 'index', 'index', 'index');
xcb = elem_centroids(mesh); uxb = unique(round(xcb(:,1)*1e9)/1e9);
tS = tic;
D = solve_static(model, [lc.F]);
t_static = toc(tS);
stat = struct([]);
for k = 1:numel(lc)
    d = D(:,k);
    [~, Rr] = solve_static(model, lc(k).F);
    rs = recover_stresses(model, d);
    s.w = d(3:5:end);
    s.twist = atand((d(wd(tipLE)) - d(wd(tipTE))) / c);
    ex = [rs.eps_top(1,:) rs.eps_bot(1,:)];
    s.ex_max = max(ex)*1e6; s.ex_min = min(ex)*1e6;
    s.tw = max(rs.tsaiwu, [], 1); s.hs = max(rs.hashin, [], 1);
    s.tw_off = interp1(uxb, arrayfun(@(u) max(s.tw(abs(xcb(:,1) - u) < 1e-9)), uxb), x_off);
    s.core = rs.core_shear; s.ilss = max(rs.ilss_index, [], 1);
    s.res = rs;
    s.eq_err = abs(sum(Rr(3:5:end)) + sum(lc(k).F(3:5:end))) / max(abs(sum(lc(k).F(3:5:end))), eps);
    stat = [stat, s]; %#ok<AGROW>
    P('  %-19s %6.2f %9.3f %9.3f %9.3f %8.3f %9.1f %9.1f %8.4f %8.4f %8.4f %8.4f\n', lc(k).name, lc(k).n, ...
      d(wd(tipLE))*1e3, d(wd(tipMC))*1e3, d(wd(tipTE))*1e3, s.twist, s.ex_max, s.ex_min, ...
      max(s.tw), max(s.hs), max(s.core), max(s.ilss));
end
P('  global equilibrium |sum(Rz) + sum(Fz)| / |sum(Fz)| (max over cases) = %.2e\n', max([stat.eq_err]));
k5 = 5;
[twmax, ewc] = max(stat(k5).tw);
[~, kpl] = max(stat(k5).res.tsaiwu(:, ewc));
hm = {'fibre tension', 'fibre compression', 'matrix tension', 'matrix compression'};
[hsmax, ehs] = max(stat(k5).hs); [~, kh] = max(stat(k5).res.hashin(:, ehs));
modeh = stat(k5).res.hashin_mode(kh, ehs);
cen = (mesh.nodes(mesh.elems(ewc,1),:) + mesh.nodes(mesh.elems(ewc,3),:))/2;
P('  ULTIMATE (LC5, n = %.2f):\n', lc(k5).n);
P('    Tsai-Wu  max %.4f -> reserve factor RF = %.2f  (ply %d, %+d deg, element at x=%.3f, y=%.3f m)\n', ...
  twmax, 1/twmax, kpl, model.sections(1).layup.theta(kpl), cen(1), cen(2));
P('    Hashin   max %.4f -> RF = %.2f, governing mode: %s (ply %d)\n', hsmax, 1/hsmax, hm{modeh}, kh);
P('    core shear index %.4f -> RF = %.2f;  interlaminar index %.4f -> RF = %.2f\n', ...
  max(stat(k5).core), 1/max(stat(k5).core), max(stat(k5).ilss), 1/max(stat(k5).ilss));
P('    Tsai-Wu at x = %.0f mm (converged design value) %.4f -> RF = %.2f\n', x_off*1e3, stat(k5).tw_off, 1/stat(k5).tw_off);
P('    -> minimum RF (corner peak incl.) = %.2f  (RF > 1 required at ultimate load)\n', ...
  min(1./[twmax hsmax max(stat(k5).core) max(stat(k5).ilss)]));

% figures
fig = new_fig(1100, 520);
plot_field(mesh, stat(4).w*1e3, sprintf('LC4 critical gust (n = %.2f): deflection w [mm]', lc(4).n), 'w [mm]', stat(4).w*20);
save_fig(fig, 'M2_static_deflection_3d');

fig = new_fig(1100, 450);
subplot(1,2,1); hold on;
cl = lines(numel(lc));
for k = 1:5
    wmc = stat(k).w(mesh.eta == 0.5 | abs(mesh.eta - 0.5) < 1e-9);
    xs = X(abs(mesh.eta - 0.5) < 1e-9);
    [xs, o] = sort(xs); plot(xs, wmc(o)*1e3, '-', 'Color', cl(k,:), 'LineWidth', 1.5);
end
grid on; xlabel('span x [m]'); ylabel('mid-chord deflection [mm]'); title('Span-wise deflection');
legend(strrep({lc(1:5).name}, '_', ' '), 'Location', 'northwest');
subplot(1,2,2); hold on;
for k = 1:5
    tw = zeros(mesh.nx+1, 1); xs = zeros(mesh.nx+1, 1);
    for i = 1:mesh.nx+1
        nle = i; nte = i + mesh.ny*(mesh.nx+1);
        tw(i) = atand((stat(k).w(nle) - stat(k).w(nte)) / c); xs(i) = X(nle);
    end
    plot(xs, tw, '-', 'Color', cl(k,:), 'LineWidth', 1.5);
end
grid on; xlabel('span x [m]'); ylabel('elastic twist [deg] (+ nose up)'); title('Span-wise twist');
save_fig(fig, 'M2_static_spanwise');

fig = new_fig(1100, 800);
subplot(2,1,1);
plot_field(mesh, stat(4).res.eps_bot(1,:)'*1e6, 'LC4: bottom-skin span-wise strain \epsilon_x [\mu\epsilon]', '\mu\epsilon');
subplot(2,1,2);
plot_field(mesh, stat(5).tw', sprintf('LC5 ultimate: max Tsai-Wu index over plies (max %.3f)', twmax), 'TW');
save_fig(fig, 'M2_static_strain_failure');

fig = new_fig(1100, 800);
subplot(2,1,1);
plot_field(mesh, stat(5).core', 'LC5 ultimate: core shear index \tau_c / S_c', 'index');
subplot(2,1,2);
plot_field(mesh, stat(5).ilss', 'LC5 ultimate: interlaminar shear index (face plies)', 'index');
save_fig(fig, 'M2_static_core_ilss');

%% C. Importance of transverse shear -------------------------------------
P('\n[C] TRANSVERSE-SHEAR CONTRIBUTION (LC4, tip mid-chord deflection; f1)\n');
lay = wing_layup(cfg);
variants = {'FSDT, k energy (used)', laminate_section(lay, 'auto');
            'FSDT, k = 5/6', laminate_section(lay, 5/6);
            'shear-rigid (~ CLT)', laminate_section(lay, 1e4)};
wref = 0;
for v = 1:size(variants,1)
    mdl = assemble_system(mesh, variants{v,2}, struct('shear', cfg.fe.shear));
    mdl = apply_bc(mdl, cfg.fe.bc);
    lcv = wing_load_cases(cfg, mdl);
    dv = solve_static(mdl, lcv(4).F); fv = solve_modal(mdl, 1);
    if v == 1, wref = dv(wd(tipMC)); f1ref = fv; end
    P('  %-24s w_tip = %8.4f mm (%6.1f %%)   f1 = %7.3f Hz (%6.1f %%)\n', variants{v,1}, ...
      dv(wd(tipMC))*1e3, 100*dv(wd(tipMC))/wref, fv, 100*fv/f1ref);
end

%% D. Modal analysis ---------------------------------------------------------
nm = cfg.an.n_modes;
tM = tic; [fr, Phi, om] = solve_modal(model, nm); t_modal = toc(tM);
P('\n[D] MODAL ANALYSIS (%d modes)\n', nm);
P('  %4s %10s %8s   %s\n', 'mode', 'f [Hz]', 'type', 'energy fractions [bend tors chord in-plane]');
labels = cell(nm,1);
for k = 1:nm
    [labels{k}, fr_k] = classify_mode(model, Phi(:,k));
    P('  %4d %10.3f %8s   [%.2f %.2f %.2f %.2f]\n', k, fr(k), labels{k}, fr_k);
end
orth = Phi' * model.M * Phi;
P('  M-orthonormality: max |Phi''M Phi - I| = %.2e\n', max(max(abs(orth - eye(nm)))));
P('  K-check: max |Phi''K Phi - diag(w^2)| / w_max^2 = %.2e\n', ...
  max(max(abs(Phi'*model.K*Phi - diag(om.^2)))) / om(end)^2);
fig = new_fig(1200, 750);
for k = 1:6
    subplot(2,3,k);
    if strcmp(labels{k}, 'in-plane')
        % in-plane (chord-wise, lead-lag) mode: w is ~0, show the chord-wise v instead
        wk = Phi(2:5:end, k); wk = wk / max(abs(wk));
        ttl = sprintf('mode %d: %.2f Hz (in-plane, v shown)', k, fr(k));
    else
        wk = Phi(3:5:end, k); wk = wk / max(abs(wk));
        ttl = sprintf('mode %d: %.2f Hz (%s)', k, fr(k), labels{k});
    end
    plot_field(mesh, wk, ttl, '', wk*0.15);
    colorbar off;
end
save_fig(fig, 'M2_mode_shapes');

%% E. Dynamic gust response -----------------------------------------------
P('\n[E] DYNAMIC 1-COS GUST RESPONSE (Newmark average acceleration)\n');
[al, be] = rayleigh_coeffs(om(cfg.an.zeta_modes(1)), om(cfg.an.zeta_modes(2)), cfg.an.zeta);
H = cfg.an.gust_H_chords * cfg.ac.cbar;
Vg = env.critical_V;
tg = 2*H/Vg;
dn = env.n_limit_pos - 1;
rz = zeros(model.ndof,1); rz(3:5:end) = 1;
F1 = (lc(1).F + cfg.g*model.M*rz);          % 1 g lift only (unit load factor)
Finc = dn*F1 - dn*cfg.g*model.M*rz;         % incremental lift - inertia relief
d1g = D(:,1);
gfun = @(t) 0.5*(1 - cos(2*pi*t/tg)) .* (t <= tg);
dt = 2.5e-4; t = 0:dt:0.6;
sensor_dofs = [wd(tipLE); wd(tipMC); wd(tipTE)];
tN = tic;
out = solve_newmark(model, Finc, gfun, t, struct('alpha', al, 'beta_r', be, 'out_dofs', sensor_dofs));
t_newmark = toc(tN);
dstat_inc = solve_static(model, Finc);
DAF = max(out.d(2,:)) / dstat_inc(wd(tipMC));
P('  gust at V = %.1f m/s, H = %.3f m, duration t_g = %.3f s (%.2f Hz), f1 = %.2f Hz, t_g/T1 = %.1f\n', ...
  Vg, H, tg, 1/tg, fr(1), tg*fr(1));
P('  Rayleigh alpha = %.4f 1/s, beta = %.3e s (zeta = %.0f %% at modes %d & %d)\n', al, be, ...
  100*cfg.an.zeta, cfg.an.zeta_modes);
P('  peak tip deflection (1g + gust) = %.3f mm, dynamic amplification factor = %.3f\n', ...
  (d1g(wd(tipMC)) + max(out.d(2,:)))*1e3, DAF);
P('  peak tip acceleration = %.2f m/s^2 (%.2f g)\n', max(abs(out.a(2,:))), max(abs(out.a(2,:)))/cfg.g);

fig = new_fig(1000, 650);
subplot(2,1,1);
plot(t, (d1g(wd(tipMC)) + out.d(2,:))*1e3, 'b-', t, (d1g(wd(tipMC)) + dstat_inc(wd(tipMC))*gfun(t))*1e3, 'r--', 'LineWidth', 1.4);
grid on; ylabel('tip deflection [mm]'); legend('dynamic (Newmark)', 'quasi-static', 'Location', 'northeast');
title(sprintf('1-cos gust at V_D: t_g = %.3f s, DAF = %.3f', tg, DAF));
subplot(2,1,2);
plot(t, out.a(2,:)/cfg.g, 'k-', 'LineWidth', 1.2); grid on;
xlabel('time [s]'); ylabel('tip acceleration [g]');
save_fig(fig, 'M2_gust_response');

% gust-duration sweep: DAF vs t_g / T1 (shock response spectrum)
ratios = [0.25 0.5 0.75 1 1.5 2 3 5 8];
DAFs = zeros(size(ratios));
for i = 1:numel(ratios)
    tgi = ratios(i)/fr(1);
    gi = @(tt) 0.5*(1 - cos(2*pi*tt/tgi)) .* (tt <= tgi);
    ti = 0:min(dt, tgi/80):(tgi + 3/fr(1));
    oi = solve_newmark(model, Finc, gi, ti, struct('alpha', al, 'beta_r', be, 'out_dofs', wd(tipMC)));
    DAFs(i) = max(oi.d) / dstat_inc(wd(tipMC));
end
P('  DAF vs t_g/T1: %s\n', sprintf('%.2f->%.3f  ', [ratios; DAFs]));
fig = new_fig(750, 480);
semilogx(ratios, DAFs, 'bo-', 'LineWidth', 1.5); hold on;
semilogx(tg*fr(1), DAF, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
grid on; xlabel('gust duration t_g / T_1'); ylabel('dynamic amplification factor');
title('Gust response spectrum of the wing (1-cos pulse, \zeta = 2%)');
legend('sweep', 'design gust', 'Location', 'northeast');
save_fig(fig, 'M2_gust_spectrum');

%% F. Damage sensitivity demo ----------------------------------------------
P('\n[F] DAMAGE SENSITIVITY (bridge to Milestone 3)\n');
dm(1) = struct('type', 'debond',      'center', [0.30 0.125], 'radius', [0.06 0.04], 'severity', 0.9, 'profile', 'uniform', 'face', 'bottom');
dm(2) = struct('type', 'skin_matrix', 'center', [0.30 0.125], 'radius', [0.06 0.04], 'severity', 0.5, 'profile', 'uniform', 'face', 'bottom');
dm(3) = struct('type', 'skin_fiber',  'center', [0.15 0.125], 'radius', [0.05 0.05], 'severity', 0.5, 'profile', 'gaussian', 'face', 'bottom');
dm(4) = struct('type', 'core_crush',  'center', [0.60 0.060], 'radius', [0.05 0.04], 'severity', 0.8, 'profile', 'uniform', 'face', 'bottom');
P('  %-12s %5s %4s  %s  %10s %8s\n', 'damage', 'D', 'nel', 'frequency shift df/f [%] modes 1..6', 'dw_tip[%]', 'min MAC');
for i = 1:numel(dm)
    [secs, sid, di] = apply_damage(cfg, mesh, dm(i));
    mdl = build_wing_model(cfg, secs, sid);
    [frd, Phd] = solve_modal(mdl, 6);
    dd = solve_static(mdl, lc(4).F);
    mac = zeros(1,6);
    for k = 1:6       % MAC on full DOF vectors (mode 4 is in-plane: w alone is noise)
        a1 = Phi(:,k); a2 = Phd(:,k);
        mac(k) = (a1'*a2)^2 / ((a1'*a1)*(a2'*a2));
    end
    P('  %-12s %5.2f %4d  [%s]  %10.4f %8.5f\n', dm(i).type, dm(i).severity, di.n_damaged, ...
      sprintf('%7.3f ', 100*(frd(:)./fr(1:6) - 1)), 100*(dd(wd(tipMC))/D(wd(tipMC),4) - 1), min(mac));
    if i == 1
        rsd = recover_stresses(mdl, dd);
        dex = (rsd.eps_bot(1,:) - stat(4).res.eps_bot(1,:))*1e6;
        fig = new_fig(1100, 800);
        subplot(2,1,1); plot_field(mesh, di.Dmax, 'Damage map: core debond D = 0.9 (bottom face)', 'D');
        subplot(2,1,2); plot_field(mesh, dex', 'LC4 strain change \Delta\epsilon_x bottom skin [\mu\epsilon] due to debond', '\mu\epsilon');
        save_fig(fig, 'M2_damage_demo');
    end
end

%% G. Cost ---------------------------------------------------------------
P('\n[G] COMPUTATIONAL COST (this machine, %s)\n', ifelse_str(exist('OCTAVE_VERSION','builtin'), 'GNU Octave', 'MATLAB'));
P('  assembly %.3f s, 6 static load cases %.3f s, %d modes %.3f s, Newmark %d steps %.2f s\n', ...
  t_asm, t_static, nm, t_modal, numel(t), t_newmark);
P('  unique element matrices computed: %d of %d elements (cache)\n', model.n_unique_elems, mesh.ne);

env = strip_handles(env);      % function handles are not portable in .mat files
save(fullfile(root, 'results', 'data', 'M2_results.mat'), 'cfg', 'conv', 'meshes', 'qinf', 'fr', ...
     'labels', 'DAF', 'DAFs', 'ratios', 'out', 'tg', 'lc', 'env', '-v7');
P('\nMilestone 2 wing analysis complete.\n');
fclose(fid);
