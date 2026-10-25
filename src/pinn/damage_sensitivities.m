function S = damage_sensitivities(spec, cfg0, ctx, sens, delta)
%DAMAGE_SENSITIVITIES  FE-based linear physics decoder: damage -> features.
%
%   S = damage_sensitivities(spec, cfg0, ctx, sens, delta)
%   For every damage TYPE t and element e, the first-order change of the
%   measurable physics features per unit damage severity D_e:
%
%       Delta x  ~=  J_t(:, e) * D_e
%
%   computed on the healthy reference wing (Milestone 4 baseline) from the
%   discretised governing equations themselves:
%
%   * natural frequencies, eigenvalue perturbation (mass unchanged):
%         d lambda_k / d D_e = phi_k' (dK_e/dD) phi_k
%         d f_k / f_k        = (d lambda_k / d D_e) / (2 lambda_k)
%     -> frequency-RATIO features  f_k/f_1 (k = 2..5)
%   * static tip-test strains, direct differentiation of  K u = F:
%         du/dD_e = -K^-1 (dK_e/dD) u
%         d eps_s = (Bm -/+ h/2 Bb) du     at every FBG, bottom/top skin
%     -> strain-SHAPE residual features (the global scale is projected out
%        exactly as in extract_features:  dx = (I - e0 e0'/e0'e0) d eps / max|e0|)
%
%   dK_e/dD is obtained by a forward difference of the ELEMENT stiffness
%   with the damage module's own degradation rules (delta = 1e-3), so the
%   decoder knows exactly what each damage type does physically.
%
%   Types: 1 debond  2 matrix-bottom  3 matrix-top  4 fibre-bottom
%          5 fibre-top  6 core crush  7 delam-bottom  8 delam-top
%   S.J (n_feat x ne x 8),  S.Jroot (n_feat x 1): uniform root-spring loss,
%   per unit (severity x loosened fraction of the root).
%   Feature order: [f2/f1 .. f5/f1 (4), eps bottom FBG (30), eps top FBG (30)]
if nargin < 5, delta = 1e-3; end
names = {spec.vars.name};
p = zeros(1, spec.k);
for j = 1:spec.k
    switch names{j}
        case {'E1','E2','G12','Gc','Ec','rho_f','rho_c','tply','tcore'}, p(j) = 1;
        case 'T', p(j) = 20;
        case 'kroot', p(j) = 10^mean([spec.vars(j).lo spec.vars(j).hi]);
    end
end
[cfg, ~] = sample_to_case(spec, p, 'healthy', cfg0);
mesh = wing_mesh(cfg.geom, cfg.fe.nx, cfg.fe.ny);
cfg.fe.root_factor = ones(numel(mesh.sets.root), 1);
model = build_wing_model(cfg);
[fr, Phi, om] = solve_modal(model, ctx.n_modes);
flex = find(arrayfun(@(k) ~strcmp(classify_mode(model, Phi(:,k)), 'in-plane'), 1:ctx.n_modes), 5);
lam = om(flex).^2; Ph = Phi(:, flex); fr = fr(flex);

% static tip test and FBG strain operator
tipy = cfg.geom.span*tand(cfg.geom.sweep_le) + 0.5*cfg.geom.c_root*cfg.geom.taper;
F2 = point_load(model.mesh, [cfg.geom.span, tipy], 3, 50);
u = solve_static(model, F2);
sec0 = model.sections(1); hs = sec0.h/2;
nf = numel(sens.fbg_elem);
E = zeros(2*nf, model.ndof);                     % eps_x at FBGs, bottom then top
for s = 1:nf
    e = sens.fbg_elem(s); d = model.edofs(e,:);
    E(s, d)    = model.rec.Bm(1,:,e) - hs*model.rec.Bb(1,:,e);
    E(nf+s, d) = model.rec.Bm(1,:,e) + hs*model.rec.Bb(1,:,e);
end
e0 = E*u;
Pproj = (eye(2*nf) - e0*e0'/(e0'*e0)) / max(abs(e0));

% element stiffness derivative per damage type (uniform mesh -> one element)
xy = model.mesh.nodes(model.mesh.elems(1,:), :);
K0 = elem_fsdt_q4(xy, sec0, cfg.fe.shear);
xc1 = mean(xy, 1);
tdef = {'debond','bottom'; 'skin_matrix','bottom'; 'skin_matrix','top'; 'skin_fiber','bottom'; ...
        'skin_fiber','top'; 'core_crush','bottom'; 'delam','bottom'; 'delam','top'};
nt = size(tdef,1); dK = zeros(20, 20, nt);
for t = 1:nt
    dm = struct('type', tdef{t,1}, 'center', xc1, 'radius', [1e-3 1e-3], 'severity', delta, ...
                'profile', 'uniform', 'face', tdef{t,2});
    [secs, sid] = apply_damage(cfg, model.mesh, dm, 1);
    Kd = elem_fsdt_q4(xy, secs(sid(1)), cfg.fe.shear);
    dK(:,:,t) = (Kd - K0) / delta;
end

ne = model.mesh.ne; nfeat = 4 + 2*nf;
S.J = zeros(nfeat, ne, nt); S.Je = zeros(2*nf, ne, nt);   % Je: raw FBG strain change per unit D
Kf = model.K(model.free, model.free);
for t = 1:nt
    Sf = zeros(5, ne); RHS = zeros(model.ndof, ne);
    for e = 1:ne
        d = model.edofs(e,:);
        Sf(:,e) = sum(Ph(d,:) .* (dK(:,:,t)*Ph(d,:)), 1)' ./ (2*lam(:));
        RHS(d, e) = -dK(:,:,t)*u(d);
    end
    du = zeros(model.ndof, ne);
    du(model.free,:) = Kf \ RHS(model.free,:);
    Jr = (fr(2:5)/fr(1)) .* (Sf(2:5,:) - Sf(1,:));
    S.Je(:,:,t) = E*du;
    S.J(:,:,t) = [Jr; Pproj*S.Je(:,:,t)];
end

% root loosening: finite difference of the whole feature vector
cfgr = cfg; cfgr.fe.root_factor = (1 - delta)*ones(numel(mesh.sets.root), 1);
mr = build_wing_model(cfgr);
[frr, Phr] = solve_modal(mr, ctx.n_modes);
ur = solve_static(mr, F2);
x0 = [fr(2:5)/fr(1); Pproj*(e0 - e0)];
xr = [frr(flex(2:5))/frr(flex(1)); Pproj*(E*ur - e0)];
S.Jroot = (xr - x0) / delta;

S.types = tdef; S.f_ref = fr; S.e0 = e0; S.elem_xy = elem_centroids(model.mesh);
S.feat_names = [arrayfun(@(k) sprintf('f%d/f1', k), 2:5, 'UniformOutput', false), ...
                arrayfun(@(k) sprintf('epsB_s%02d', k), 1:nf, 'UniformOutput', false), ...
                arrayfun(@(k) sprintf('epsT_s%02d', k), 1:nf, 'UniformOutput', false)];
S.delta = delta;
end
