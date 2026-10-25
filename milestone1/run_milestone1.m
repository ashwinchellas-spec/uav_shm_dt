%RUN_MILESTONE1  Step 1 - Problem definition & wing modelling.
%
%   Produces every number and figure quoted in reports/Milestone1_Report.md:
%     1. configuration summary
%     2. equivalent sandwich depth of the NACA 0012 wing
%     3. laminate (CLT + FSDT) section properties and shear correction
%     4. flight envelope (V-n) and design load factors
%     5. aerodynamic load model and load-case resultants
%     6. finite-element mesh, boundary conditions and model size
%   Outputs: results/figures/M1_*.png, results/data/M1_summary.txt,
%            results/data/M1_baseline.mat
%
%   Run from the project root:  >> run_milestone1

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
cfg  = wing_baseline_config();
fid  = fopen(fullfile(root, 'results', 'data', 'M1_summary.txt'), 'w');
P = @(varargin) logprintf(fid, varargin{:});

P('=====================================================================\n');
P(' MILESTONE 1 - Problem definition: %s\n', cfg.name);
P('=====================================================================\n');

%% 1. Geometry ---------------------------------------------------------
g = cfg.geom;
P('\n[1] GEOMETRY\n');
P('  semi-span b = %.3f m, root chord = %.3f m, taper = %.2f, LE sweep = %.1f deg\n', ...
  g.span, g.c_root, g.taper, g.sweep_le);
P('  airfoil %s (t/c = %.2f), full wing area S = %.3f m^2, AR = %.2f\n', ...
  g.airfoil, g.tc, cfg.ac.S, cfg.ac.b_full^2/cfg.ac.S);

%% 2. Equivalent sandwich depth ----------------------------------------
% Skins of thickness t following the airfoil contour have (per unit
% width of skin, about the chord line)  I = 2 t int_0^c y_t(x)^2 dx.
% A flat sandwich with face-centroid distance d has I = 2 t c (d/2)^2.
% Equal I  ->  d_eq = 2 sqrt( mean(y_t^2) ).
s = linspace(0, 1, 20001);
yt = naca4_thickness(s, g.tc);
d_eq = 2*sqrt(trapz(s, yt.^2)) * g.c_root;
t_face = cfg.lam.t_ply * numel(cfg.lam.face_angles);
d_model = cfg.lam.t_core + t_face;
P('\n[2] EQUIVALENT SANDWICH DEPTH (NACA 0012)\n');
P('  max thickness t_max = %.2f mm, mean thickness = %.2f mm\n', ...
  g.tc*g.c_root*1e3, 2*trapz(s, yt)*g.c_root*1e3);
P('  bending-equivalent face distance d_eq = %.2f mm\n', d_eq*1e3);
P('  model: core %.2f mm + face %.3f mm -> d = %.2f mm (%.1f %% of d_eq)\n', ...
  cfg.lam.t_core*1e3, t_face*1e3, d_model*1e3, 100*d_model/d_eq);

fig = new_fig(900, 380);
c = g.c_root*1e3; hold on;
plot(s*c, yt*c, 'b-', s*c, -yt*c, 'b-', 'LineWidth', 1.5);
plot([0 c], [d_model d_model]/2*1e3, 'r--', [0 c], -[d_model d_model]/2*1e3, 'r--', 'LineWidth', 1.5);
axis equal; grid on; xlabel('chord y [mm]'); ylabel('z [mm]');
title('NACA 0012 skin contour (blue) vs. equivalent sandwich face planes (red)');
legend('airfoil upper skin', 'airfoil lower skin', 'equivalent faces', 'Location', 'southoutside');
save_fig(fig, 'M1_equivalent_sandwich');

%% 3. Laminate ---------------------------------------------------------
layup = wing_layup(cfg);
sec   = laminate_section(layup, cfg.lam.kshear);
sec56 = laminate_section(layup, 5/6);
P('\n[3] LAMINATE (bottom -> top)\n  ');
for k = 1:numel(layup.theta)
    if strcmp(layup.mat(k).type, 'core'), P('[CORE %.2fmm] ', layup.t(k)*1e3);
    else, P('%+d ', layup.theta(k)); end
end
P('\n  total thickness h = %.3f mm, areal mass I0 = %.3f kg/m^2\n', sec.h*1e3, sec.I0);
semi_mass = sec.I0 * cfg.ac.S/2;
P('  semi-wing structural mass = %.3f kg (%.1f %% of MTOW for both wings)\n', ...
  semi_mass, 200*semi_mass/cfg.ac.mass);
P('  A [MN/m]:\n'); P('    %10.4f %10.4f %10.4f\n', sec.A'/1e6);
P('  B [N]   : max |B_ij| = %.3e  (symmetric laminate -> 0)\n', max(abs(sec.B(:))));
P('  D [N m]:\n');  P('    %10.4f %10.4f %10.4f\n', sec.D');
P('  As (k energy) [MN/m]:  [%.4f %.4f; %.4f %.4f]\n', sec.As'/1e6);
P('  shear correction (energy method): kx = %.4f, ky = %.4f   (homogeneous: 0.8333)\n', ...
  sec.k_energy);

% sandwich-theory cross-check (Allen 1969)
fl.mat = repmat(material_library(cfg.lam.face_mat), 1, numel(cfg.lam.face_angles));
fl.theta = cfg.lam.face_angles; fl.t = cfg.lam.t_ply*ones(1, numel(fl.theta));
fs = laminate_section(fl, 5/6);
af = inv(fs.A);
Ef = [1/(t_face*af(1,1)), 1/(t_face*af(2,2))];
D_sw  = Ef * t_face * d_model^2 / 2;        % thin-face sandwich bending
cm = material_library(cfg.lam.core_mat);
S_sw  = cm.G13 * d_model^2 / cfg.lam.t_core; % sandwich shear stiffness
P('  face membrane moduli Ex = %.2f GPa, Ey = %.2f GPa (quasi-isotropic)\n', Ef/1e9);
P('  sandwich theory D11 = %.2f, D22 = %.2f N m  vs CLT %.2f, %.2f N m\n', D_sw, sec.D(1,1), sec.D(2,2));
P('  sandwich theory shear stiffness  G_c d^2/t_c = %.4f MN/m\n', S_sw/1e6);
P('     FSDT, k from energy method    As11        = %.4f MN/m  (ratio %.3f)\n', sec.As(1,1)/1e6, sec.As(1,1)/S_sw);
P('     FSDT, k = 5/6                 As11        = %.4f MN/m  (ratio %.1f  <-- wrong)\n', sec56.As(1,1)/1e6, sec56.As(1,1)/S_sw);

% layup stack figure (schematic, not to scale: plies drawn 1 unit, core 5 units)
fig = new_fig(800, 520); hold on;
cols = [0.20 0.40 0.80; 0.90 0.50 0.10; 0.30 0.70 0.30; 0.60 0.20 0.60];
angs = [0 45 -45 90];
ybot = 0;
for k = 1:numel(layup.theta)
    if strcmp(layup.mat(k).type, 'core')
        hd = 5; fc = [1 0.9 0.6];
        lab = sprintf('Rohacell 51 WF foam core, t = %.2f mm', layup.t(k)*1e3);
    else
        hd = 1; fc = cols(find(angs == layup.theta(k), 1), :);
        lab = sprintf('ply %d: T300/5208, %+d deg, t = %.3f mm', k, layup.theta(k), layup.t(k)*1e3);
    end
    patch([0 1 1 0], [ybot ybot ybot+hd ybot+hd], fc, 'EdgeColor', 'k');
    text(1.05, ybot + hd/2, lab, 'FontSize', 10);
    ybot = ybot + hd;
end
xlim([0 2.3]); ylim([-0.5 ybot+0.5]); axis off;
title(sprintf('Sandwich stacking [0/45/-45/90/C]_s  (schematic, not to scale; h = %.2f mm)', sec.h*1e3));
text(0.02, -0.3, 'bottom (z = -h/2)', 'FontSize', 9); text(0.02, ybot+0.3, 'top (z = +h/2)', 'FontSize', 9);
save_fig(fig, 'M1_layup');

% through-thickness transverse shear stress profile
fig = new_fig(900, 450);
subplot(1,2,1);
plot(sec.tprof.fx, sec.tprof.z*1e3, 'b-', 'LineWidth', 1.5); hold on;
plot(ones(1,2)/sec.h, [sec.z(1) sec.z(end)]*1e3, 'r--', 'LineWidth', 1.2);
grid on; xlabel('\tau_{xz}(z) / Q_x  [1/m]'); ylabel('z [mm]');
title('Sandwich: equilibrium profile vs FSDT');
legend('equilibrium (used)', 'FSDT constant', 'Location', 'southoutside');
subplot(1,2,2);
m1 = material_library('T300_5208'); hl.mat = repmat(m1,1,8);
hl.theta = [0 45 -45 90 90 -45 45 0]; hl.t = 0.125e-3*ones(1,8);
hs = laminate_section(hl, 'auto');
plot(hs.tprof.fx, hs.tprof.z*1e3, 'b-', 'LineWidth', 1.5); hold on;
plot(ones(1,2)/hs.h, [hs.z(1) hs.z(end)]*1e3, 'r--', 'LineWidth', 1.2);
grid on; xlabel('\tau_{xz}(z) / Q_x  [1/m]'); ylabel('z [mm]');
title(sprintf('Solid [0/45/-45/90]_s: k_x = %.3f', hs.k_energy(1)));
save_fig(fig, 'M1_shear_profiles');
P('  solid 8-ply [0/45/-45/90]s laminate for comparison: kx = %.4f, ky = %.4f\n', hs.k_energy);

%% 4. Flight envelope ----------------------------------------------------
env = flight_envelope(cfg.ac, cfg.g);
P('\n[4] FLIGHT ENVELOPE (V-n)\n');
P('  W = %.2f N, W/S = %.1f N/m^2, V_S1 = %.2f m/s, V_A = %.2f m/s, V_C = %.1f, V_D = %.1f m/s\n', ...
  cfg.ac.mass*cfg.g, env.gust.WS, env.Vs1, env.Va, cfg.ac.Vc, cfg.ac.Vd);
P('  lift slope a = %.3f /rad, mass ratio mu = %.2f, gust alleviation Kg = %.3f\n', ...
  env.gust.a, env.gust.mu, env.gust.Kg);
P('  gust dn at V_C (%.2f m/s) = %.2f, at V_D (%.2f m/s) = %.2f\n', ...
  cfg.ac.Ude_c, env.dn_c, cfg.ac.Ude_d, env.dn_d);
for k = 1:numel(env.points)
    P('    %-26s V = %5.1f m/s   n = %+6.2f\n', env.points(k).name, env.points(k).V, env.points(k).n);
end
P('  >> design limit load factors: n+ = %.2f (%s),  n- = %.2f (%s)\n', ...
  env.n_limit_pos, env.critical, env.n_limit_neg, env.critical_neg);

fig = new_fig(950, 620); hold on;
% manoeuvre envelope as one closed polyline
Va = env.Va; Vd = cfg.ac.Vd; Vc = cfg.ac.Vc;
V1 = linspace(0, Va, 100);
Vneg = sqrt(2*cfg.ac.mass*cfg.g*cfg.ac.n_neg/(cfg.ac.rho*cfg.ac.S*cfg.ac.CLmin));  % negative corner
V2 = linspace(0, Vneg, 100);
envV = [V1, Vd, Vd, Vc, fliplr([V2 Vneg])];
envN = [env.nstall_p(V1), cfg.ac.n_pos, 0, cfg.ac.n_neg, fliplr([env.nstall_n(V2) cfg.ac.n_neg])];
envV = [envV(1:numel(V1)), Va, envV(numel(V1)+1:end)]; envN = [envN(1:numel(V1)), cfg.ac.n_pos, envN(numel(V1)+1:end)];
h1 = plot(envV, envN, 'b-', 'LineWidth', 2.2);
Vs = linspace(0, Vd, 200);
h2 = plot(Vs, env.nstall_p(Vs), 'b:', 'LineWidth', 1);
plot(Vs, env.nstall_n(Vs), 'b:', 'LineWidth', 1);
h3 = plot([0 Vc], [1 1+env.dn_c], 'r--', 'LineWidth', 1.3);
plot([0 Vc], [1 1-env.dn_c], 'r--', 'LineWidth', 1.3);
h4 = plot([0 Vd], [1 1+env.dn_d], 'm--', 'LineWidth', 1.3);
plot([0 Vd], [1 1-env.dn_d], 'm--', 'LineWidth', 1.3);
pv = [env.points.V]; pn = [env.points.n];
h5 = plot(pv, pn, 'ks', 'MarkerSize', 7, 'MarkerFaceColor', [0.6 0.6 0.6]);
h6 = plot(env.critical_V, env.n_limit_pos, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
grid on; xlabel('equivalent airspeed V [m/s]'); ylabel('load factor n');
title(sprintf('V-n diagram, %g kg UAV: design limit n+ = %.2f, n- = %.2f', cfg.ac.mass, env.n_limit_pos, env.n_limit_neg));
legend([h1 h2 h3 h4 h5 h6], {'manoeuvre envelope', 'stall boundary (C_{Lmax}, C_{Lmin})', ...
    sprintf('gust U_{de} = %.2f m/s (V_C)', cfg.ac.Ude_c), sprintf('gust U_{de} = %.2f m/s (V_D)', cfg.ac.Ude_d), ...
    'design points (stall-capped)', 'critical point'}, 'Location', 'northwest');
ylim([-4 9]); xlim([0 Vd*1.05]);
save_fig(fig, 'M1_Vn_diagram');

%% 5. Loads --------------------------------------------------------------
model = build_wing_model(cfg);
[lc, env] = wing_load_cases(cfg, model);
mesh = model.mesh;
xw = mesh.nodes(:,1); yw = mesh.nodes(:,2);
y_ea = 0.5*g.c_root;      % reference axis for torque: mid-chord
P('\n[5] LOAD CASES (root reactions, from the consistent nodal loads)\n');
P('  %-20s %7s %10s %12s %12s\n', 'case', 'n', 'Fz [N]', 'Mx_root[Nm]', 'T(mid-c)[Nm]');
for k = 1:numel(lc)
    Fw = lc(k).F(3:5:end);
    lc(k).Fz = sum(Fw);
    lc(k).Mroot = sum(Fw .* xw);
    lc(k).T = sum(Fw .* (y_ea - yw));
    P('  %-20s %7.2f %10.2f %12.3f %12.3f\n', lc(k).name, lc(k).n, lc(k).Fz, lc(k).Mroot, lc(k).T);
end
P('  (torque positive nose-up; lift acts at quarter chord, ahead of mid-chord)\n');

Wt = cfg.ac.mass*cfg.g;
[~, Le] = lift_distribution(g, Wt/2, 'elliptic');
[~, Lp] = lift_distribution(g, Wt/2, 'planform');
[~, Ls] = lift_distribution(g, Wt/2, 'schrenk');
xs = linspace(0, g.span, 300);
fig = new_fig(1000, 420);
subplot(1,2,1);
plot(xs, Le(xs), 'r--', xs, Lp(xs), 'g--', xs, Ls(xs), 'b-', 'LineWidth', 1.6);
grid on; xlabel('span x [m]'); ylabel('L''(x) [N/m] at n = 1');
legend('elliptic', 'planform', 'Schrenk (used)', 'Location', 'southwest');
title('Span-wise lift distribution');
subplot(1,2,2);
ss = linspace(0.005, 1, 300);
plot(ss, (2/pi)*sqrt((1-ss)./ss), 'b-', 'LineWidth', 1.6); hold on;
plot([0.25 0.25], [0 4], 'k--');
grid on; xlabel('s = y / c'); ylabel('f(s)'); ylim([0 4]);
title('Chord-wise loading (thin airfoil); x_{cp} = c/4');
save_fig(fig, 'M1_lift_distribution');

%% 6. Mesh & model size --------------------------------------------------
P('\n[6] FINITE-ELEMENT MODEL\n');
P('  Q4 FSDT elements: %d x %d = %d, nodes = %d, DOFs = %d (free %d)\n', ...
  mesh.nx, mesh.ny, mesh.ne, mesh.nn, model.ndof, numel(model.free));
P('  element size %.1f x %.1f mm, element a/h = %.2f, shear = %s, BC = %s\n', ...
  1e3*g.span/mesh.nx, 1e3*g.c_root/mesh.ny, (g.span/mesh.nx)/sec.h, cfg.fe.shear, cfg.fe.bc);
rz = zeros(model.ndof, 1); rz(3:5:end) = 1;
P('  semi-wing mass from consistent M (r_z'' M r_z): %.4f kg\n', rz' * model.M * rz);

fig = new_fig(1000, 380); hold on;
for e = 1:mesh.ne
    xy = mesh.nodes(mesh.elems(e,:), :);
    patch(xy(:,1), xy(:,2), [0.85 0.9 1], 'EdgeColor', [0.3 0.3 0.6]);
end
plot(xw(mesh.sets.root), yw(mesh.sets.root), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot([0 g.span], [0.25 0.25]*g.c_root, 'r--', 'LineWidth', 1.5);
axis equal; axis([-0.05 g.span+0.05 -0.05 g.c_root+0.05]);
xlabel('span x [m]'); ylabel('chord y [m]  (LE at y = 0)');
title(sprintf('FE mesh %dx%d Q4-FSDT-HR, clamped root (squares), quarter-chord line (dashed)', mesh.nx, mesh.ny));
set(gca, 'YDir', 'reverse');
save_fig(fig, 'M1_mesh');

env = strip_handles(env);      % function handles are not portable in .mat files
save(fullfile(root, 'results', 'data', 'M1_baseline.mat'), 'cfg', 'sec', 'env', 'lc', '-v7');
P('\nMilestone 1 complete. Figures in results/figures/M1_*.png\n');
fclose(fid);
