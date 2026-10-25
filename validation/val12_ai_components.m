function R = val12_ai_components()
%VAL12_AI_COMPONENTS  Verification of the Milestone 5 AI building blocks.
%   1. PI-Net loss gradient (all 7 terms incl. the physics decoder and the
%      evidence-map gain) vs central finite differences
%   2. surrogate MLP gradient vs finite differences
%   3. reduced-order beam model (hybrid surrogate prior) vs the FE wing:
%      first bending frequency and 1 g tip deflection, baseline design
rand('state', 1); randn('state', 1); %#ok<RAND>
n = 12; ne = 9; nf = 4;
[th, ar] = pinet_init(6, [7 5], 7, ne, 3); th = th + 0.3*randn(size(th));
X = randn(n, 6);
T = struct('y', [1 2 3 4 5 6 7 2 3 4 6 1]', 'lab', logical([1 1 1 1 1 1 1 0 1 1 1 0])', 'face', nan(n,1), ...
    'qtrue', rand(n, ne), 'loc', rand(n, 2), 'ext', rand(n,1), 'rf', randn(n,1), 'dXn', randn(n, nf), 'E', randn(n, ne));
T.face([3 4 9 10 11]) = [0 1 1 0 1]; T.qtrue = T.qtrue ./ sum(T.qtrue, 2); T.loc([1 7 12],:) = nan;
ph = struct('Jn', randn(nf, ne, 8), 'Jrootn', randn(nf, 1), 's_scale', 2, 'XY', rand(ne, 2), 'ell', 0.3);
w = struct('cls',1,'face',0.7,'map',0.5,'loc',0.8,'ext',0.6,'rf',0.4,'phys',1.3,'l2',1e-3);
g1 = grad_check(@(t) pinet_lossgrad(t, ar, X, T, ph, w), th, 80);
[tm, am] = mlp_init(5, [6 4], 3, 1);
Xm = randn(9,5); Ym = randn(9,3);
g2 = grad_check(@(t) mlp_lossgrad(t, am, Xm, Ym, 1e-3), tm + 0.2*randn(size(tm)), 60);
% beam prior vs FE (baseline, elastic root k* = 1000)
cfg = wing_baseline_config(); cfg.fe.bc = 'elastic_root'; cfg.fe.kroot = 1e3;
model = build_wing_model(cfg);
[fr, Phi] = solve_modal(model, 6);
kb = find(arrayfun(@(k) strcmp(classify_mode(model, Phi(:,k)), '1B'), 1:6), 1);
d = solve_static(model, wing_net_1g(cfg, model));
tip = model.mesh.sets.tip; [~, j] = min(abs(model.mesh.eta(tip) - 0.5));
pr = beam_prior(cfg, 7.67, 0.05);
R = [vres('V12 AI components', 'PI-Net gradient, max rel. error vs finite differences', g1, 0, 1e-5), ...
     vres('V12 AI components', 'surrogate MLP gradient, max rel. error', g2, 0, 1e-5), ...
     vres('V12 AI components', 'beam/CLT prior f_1B vs FE (baseline)', pr.f1B, fr(kb), 0.05), ...
     vres('V12 AI components', 'beam/CLT prior 1 g tip deflection vs FE', pr.wtip, d(5*(tip(j)-1)+3), 0.05)];
R(1).err = g1; R(1).pass = g1 < R(1).tol; R(2).err = g2; R(2).pass = g2 < R(2).tol;
end

function F1 = wing_net_1g(cfg, model)
% net 1 g load exactly as in run_sample (Schrenk lift W/2 minus wing weight)
p1 = lift_distribution(cfg.geom, 1, 'schrenk');
F1u = pressure_load(model.mesh, p1, 4, true); F1u = F1u / sum(F1u(3:5:end));
rz = zeros(model.ndof,1); rz(3:5:end) = 1;
F1 = cfg.ac.mass*cfg.g/2*F1u - cfg.g*model.M*rz;
end
