function o = run_sample(cfg, damages, ctx)
%RUN_SAMPLE  One high-fidelity simulation of the DOE -> fixed output schema.
%
%   o = run_sample(cfg, damages, ctx)      (ctx from doe_context.m)
%
%   Analyses per sample
%     1. model build (perturbed laminate, elastic root, damaged sections)
%     2. modal analysis (n_modes): frequencies, mode types, mass-normalised
%        w-mode shapes at all nodes, bottom-skin modal strains per element
%     3. two static unit load cases (linear -> any load factor by scaling)
%          L1: net 1 g flight load (Schrenk lift - wing weight)
%          L2: 50 N ground-test tip load at mid-chord
%        deflection field, top/bottom surface strains, tip deflection, twist
%     4. residual strength at ultimate load (n_ult x L1): Tsai-Wu, Hashin,
%        core-shear and interlaminar reserve factors (x >= x_off)
%     5. gust response by exact modal superposition (+ static correction):
%        design 1-cos gust and a sharp gust (t_g = T1): DAF, peaks
%     6. quality checks: equilibrium and eigen-residual
%
%   All large arrays are stored in single precision.
t0 = tic;
o = empty_output(ctx);
try
    mesh0 = wing_mesh(cfg.geom, cfg.fe.nx, cfg.fe.ny);
    rootfac = ones(numel(mesh0.sets.root), 1);
    if isempty(damages)
        sections = laminate_section(wing_layup(cfg), cfg.lam.kshear);
        sid = ones(mesh0.ne, 1); Dmap = zeros(mesh0.ne, 1);
    else
        [sections, sid, di] = apply_damage(cfg, mesh0, damages);
        Dmap = di.Dmax;
        for j = 1:numel(damages)
            if strcmp(damages(j).type, 'root_loosening')
                rootfac = rootfac .* root_loosening_factor(mesh0, damages(j));
            end
        end
    end
    cfg.fe.root_factor = rootfac;
    model = build_wing_model(cfg, sections, sid);
    mesh = model.mesh; nm = ctx.n_modes;

    % --- modal --------------------------------------------------------
    [fr, Phi, om] = solve_modal(model, nm);
    code = zeros(1, nm); ord = zeros(1, nm);
    for k = 1:nm
        lab = classify_mode(model, Phi(:,k));
        if strcmp(lab, 'in-plane'), code(k) = 4; ord(k) = 1;
        else
            code(k) = find('BTC' == lab(end)); ord(k) = str2double(lab(1:end-1));
        end
    end
    Em = zeros(mesh.ne, 3, nm);
    hb = -sections(1).h/2;
    for e = 1:mesh.ne
        de = Phi(model.edofs(e,:), :);
        Em(e,:,:) = reshape(model.rec.Bm(:,:,e)*de + hb*model.rec.Bb(:,:,e)*de, 1, 3, nm);
    end
    fd = model.free;                      % residual on free DOFs (fixed rows hold reactions)
    KP = model.K(fd,:)*Phi;
    res_eig = max(sqrt(sum((KP - model.M(fd,:)*Phi*diag(om.^2)).^2, 1)) ./ sqrt(sum(KP.^2, 1)));

    % --- static unit loads ------------------------------------------------
    rz = zeros(model.ndof,1); rz(3:5:end) = 1;
    if isfield(ctx, 'F1unit') && numel(ctx.F1unit) == model.ndof && ctx.same_geom
        F1u = ctx.F1unit;
    else
        p1 = lift_distribution(cfg.geom, 1, 'schrenk');
        F1u = pressure_load(mesh, p1, 4, true); F1u = F1u / sum(F1u(3:5:end));
    end
    W = cfg.ac.mass*cfg.g;
    F1 = W/2*F1u - cfg.g*model.M*rz;
    tipy = cfg.geom.span*tand(cfg.geom.sweep_le) + 0.5*cfg.geom.c_root*cfg.geom.taper;
    F2 = point_load(mesh, [cfg.geom.span, tipy], 3, 50);
    FF = [F1 F2];
    [d, Rr] = solve_static(model, FF);
    eqe = max(abs(sum(Rr(3:5:end,:), 1) + sum(FF(3:5:end,:), 1)) ./ abs(sum(FF(3:5:end,:), 1)));
    tip = mesh.sets.tip; [~, oo] = sort(mesh.eta(tip)); tip = tip(oo);
    [~, jm] = min(abs(mesh.eta(tip) - 0.5));
    tipn = [tip(1) tip(jm) tip(end)];
    ctip = cfg.geom.c_root*cfg.geom.taper;
    epsb = zeros(mesh.ne, 3, 2); epst = epsb;
    hs = sections(1).h/2;
    for e = 1:mesh.ne
        de = d(model.edofs(e,:), :);
        e0 = model.rec.Bm(:,:,e)*de; ka = model.rec.Bb(:,:,e)*de;
        epsb(e,:,:) = reshape(e0 - hs*ka, 1, 3, 2);
        epst(e,:,:) = reshape(e0 + hs*ka, 1, 3, 2);
    end

    % --- residual strength at ultimate load ------------------------------
    rs = recover_stresses(model, ctx.n_ult*d(:,1));
    xc = elem_centroids(mesh); off = xc(:,1) >= ctx.x_off;
    tw = max(rs.tsaiwu, [], 1); hsn = max(rs.hashin, [], 1);
    il = max(rs.ilss_index, [], 1); cs = rs.core_shear;
    idx = max([tw; hsn; il; cs], [], 1);
    [imax, ecrit] = max(idx .* off');
    [~, kh] = max(rs.hashin(:, ecrit));

    % --- gust response (exact modal superposition) -----------------------
    zeta = ctx.alpha./(2*om) + ctx.beta*om/2;
    gl = Phi' * F1;                               % modal force, 1 g net load
    wd = 5*(tipn(2)-1) + 3;
    rootcol = find(abs(xc(:,1) - min(xc(:,1))) < 1e-9);       % first element column
    [~, j2] = min(abs(xc(rootcol,2) - mean(xc(rootcol,2))));
    er = rootcol(j2);                                          % root element, mid-chord
    Phi_s = [Phi(wd,:); squeeze(Em(er,1,:))'];
    st = [d(wd,1); epsb(er,1,1)];
    resid = st - Phi_s*(gl./om.^2);
    dts = ctx.dt;
    tg = ctx.tg; tt = 0:dts:(tg + 3/fr(1)); gf = 0.5*(1-cos(2*pi*tt/tg)).*(tt <= tg);
    g1 = modal_response(om, zeta, Phi_s, gl, gf, dts, resid);
    T1 = 1/fr(1); dts2 = T1/100; tt2 = 0:dts2:4*T1; gf2 = 0.5*(1-cos(2*pi*tt2/T1)).*(tt2 <= T1);
    g2 = modal_response(om, zeta, Phi_s, gl, gf2, dts2, resid);

    % --- pack ---------------------------------------------------------------
    o.ok = true;
    o.freq = fr(:)'; o.mode_type = code; o.mode_order = ord; o.zeta = zeta(:)';
    o.mass = rz'*model.M*rz;
    o.Wmodes = single(Phi(3:5:end, :));
    o.Emodes = single(Em);
    o.gam_L1 = gl(:)'; o.gam_L2 = (Phi'*F2)';
    o.w_static = single(d(3:5:end, :));
    o.eps_bot = single(epsb); o.eps_top = single(epst);
    o.tipw = d(5*(tipn-1)+3, :);
    o.twist = atand((d(5*(tipn(1)-1)+3,:) - d(5*(tipn(3)-1)+3,:)) / ctip);
    o.RF = 1 ./ [max(tw(off)) max(hsn(off)) max(cs(off)) max(il(off)) imax max(idx)];
    o.crit_elem = ecrit; o.crit_hashin_mode = rs.hashin_mode(kh, ecrit);
    o.gust = [max(g1.d(1,:))/st(1), max(g1.d(1,:)), max(abs(g1.a(1,:))), max(g1.d(2,:)), ...
              max(g2.d(1,:))/st(1)];
    o.Dmap = single(Dmap(:)');
    o.root_factor = rootfac(:)';
    o.qc = [eqe, res_eig];
    % never store silent NaN/Inf: a sample with non-finite results is a failure
    chk = [o.freq, o.mass, o.tipw(:)', o.twist, o.RF(1:5), o.gust, o.qc];
    if any(~isfinite(chk)) || any(~isfinite(o.Wmodes(:))) || any(~isfinite(o.Emodes(:)))
        o.ok = false; o.msg = 'non-finite results (e.g. eigensolver did not converge)';
    end
catch err
    o.ok = false; o.msg = err.message;
end
o.time = toc(t0);
end

function o = empty_output(ctx)
nm = ctx.n_modes; ne = ctx.ne; nn = ctx.nn;
o = struct('ok', false, 'msg', '', 'time', 0, 'freq', nan(1,nm), 'mode_type', zeros(1,nm), ...
    'mode_order', zeros(1,nm), 'zeta', nan(1,nm), 'mass', nan, ...
    'Wmodes', single(nan(nn,nm)), 'Emodes', single(nan(ne,3,nm)), ...
    'gam_L1', nan(1,nm), 'gam_L2', nan(1,nm), 'w_static', single(nan(nn,2)), ...
    'eps_bot', single(nan(ne,3,2)), 'eps_top', single(nan(ne,3,2)), ...
    'tipw', nan(3,2), 'twist', nan(1,2), 'RF', nan(1,6), 'crit_elem', 0, ...
    'crit_hashin_mode', 0, 'gust', nan(1,5), 'Dmap', single(zeros(1,ne)), ...
    'root_factor', ones(1, ctx.ny+1), 'qc', nan(1,2));
end
