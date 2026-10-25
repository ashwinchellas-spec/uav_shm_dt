%RUN_MILESTONE5  Step 6: physics-informed AI model development.
%
%   PART A - damage diagnosis network (SHM tier, 8400 noisy records)
%     1. physics decoder: FE sensitivities of the measurable features to
%        element damage (eigenvalue perturbation + adjoint static analysis),
%        validated against noise-free damaged/healthy twins
%     2. exact gradients of every network (finite-difference check)
%     3. PI-Net: multi-task network (class, damaged face, damage map,
%        location, extent, residual strength) trained with a PHYSICS
%        loss (its diagnosis must reproduce the measured feature changes
%        through the FE decoder) and a physics evidence-map prior
%     4. comparison with the Milestone 4 baselines on the untouched test set
%     5. ablations: evidence prior (none / FBG / dense), physics-loss weight
%        x labelled fraction (label efficiency, 2 seeds)
%   PART B - hybrid physics-AI design surrogate (design tier, 400 wings)
%     6. reduced-order beam/CLT model + neural correction of its model-form
%        error vs a pure neural network: learning curves, extrapolation
%   All networks, losses and gradients are written from scratch (no toolbox).
%
%   >> setup_paths; run_milestone5            (~15 min Octave, less in MATLAB)

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
QUICK  = false;                 % true: fewer epochs / seeds (smoke test only)
LAMBDA = 0.3;                   % physics-loss weight of the final model
HIDDEN = [64 64];
EPOCHS = 300; if QUICK, EPOCHS = 30; end
outdir = fullfile(root, 'results', 'models');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(fullfile(root, 'results', 'data', 'M5_summary.txt'), 'w');
P_ = @(varargin) logprintf(fid, varargin{:});
P_('=====================================================================\n');
P_(' MILESTONE 5 - Physics-informed AI model development\n');
P_('=====================================================================\n');
t_all = tic;

F = load(fullfile(root, 'results', 'features', 'shm_features.mat'));
db = load(fullfile(root, 'results', 'doe', 'shm_database.mat'), 'spec', 'root_factor', 'P', 'Dmap', 'cls', 'classes');
cfg0 = wing_baseline_config();
ctx = doe_context(db.spec, cfg0);
y = F.y_class; sid = F.sample_id; nc = numel(F.class_names);
cn = strrep(F.class_names, '_', ' ');
mesh = wing_mesh(cfg0.geom, cfg0.fe.nx, cfg0.fe.ny);
sub = @(T, m) structfun(@(v) v(m,:), T, 'UniformOutput', false);

%% 1. Physics decoder ---------------------------------------------------------
P_('\n[1] PHYSICS DECODER (FE sensitivities, governing equations)\n');
tic; S = damage_sensitivities(db.spec, cfg0, ctx, F.sensors);
P_('  J: %d features x %d elements x %d damage types + root loosening (%.1f s)\n', size(S.J,1), size(S.J,2), size(S.J,3), toc);
P_('  features: f2/f1..f5/f1 (eigenvalue perturbation), 60 FBG strain-shape residuals (adjoint statics)\n');
tw = load(fullfile(root, 'results', 'doe', 'twins.mat'));
fe = F.sensors.fbg_elem;
P_('  validation against noise-free nonlinear twins (same wing with / without damage):\n');
P_('  %-14s %4s %18s %16s %16s\n', 'class', 'n', 'cosine(pred,true)', '|pred|/|true|', 'true |d eps|');
fig = new_fig(1100, 300); dec_tab = [];
for c = 2:6
    q = find(db.cls(tw.idx) == c); cc = []; ra = []; sg = []; pt = []; tt = [];
    for j = q'
        i = tw.idx(j); Dd = double(db.Dmap(i,:))'; top = db.P(i,18) > 0.5;
        switch c, case 2, t = 1; case 3, t = 2 + top; case 4, t = 4 + top; case 5, t = 6; case 6, t = 7 + top; end
        true_de = tw.eps_d(fe,j) - tw.eps_h(fe,j);
        pred_de = S.Je(1:30,:,t) * Dd;
        if norm(true_de) > 0
            cc(end+1) = (true_de'*pred_de)/(norm(true_de)*norm(pred_de) + eps); %#ok<AGROW>
            ra(end+1) = norm(pred_de)/norm(true_de); sg(end+1) = 1e6*norm(true_de); %#ok<AGROW>
            pt = [pt; 1e6*pred_de]; tt = [tt; 1e6*true_de]; %#ok<AGROW>
        end
    end
    if median(sg) < 0.01
        P_('  %-14s %4d   no static strain change (max %.3f ue): invisible to the FBG tip test\n', cn{c}, numel(q), max([sg 0]));
        continue
    end
    P_('  %-14s %4d %11.3f (p10 %.3f) %13.2f %13.2f ue\n', cn{c}, numel(cc), median(cc), prctile_simple(cc, 10), median(ra), median(sg));
    dec_tab(end+1,:) = [c numel(cc) median(cc) prctile_simple(cc,10) median(ra) median(sg)]; %#ok<AGROW>
    subplot(1,4,c-1); plot(tt, pt, '.', 'MarkerSize', 6); hold on;
    lim = max(abs([tt; pt])) + eps; plot([-lim lim], [-lim lim], 'k--');
    axis([-lim lim -lim lim]); axis square; grid on;
    xlabel('true \Delta\epsilon (\mu\epsilon)'); ylabel('decoder (\mu\epsilon)'); title(cn{c});
end
save_fig(fig, 'M5_decoder_validation');
P_('  -> the decoder reproduces the SHAPE of the strain change (cosine ~1); linearisation\n');
P_('     under-predicts large damage by 5-25 %%.  The true signal is only ~1-3 ue (FBG noise 1 ue + 0.5 %% gain + 1 %% load error)\n');

%% 2. Verification of the AI building blocks -------------------------------------
P_('\n[2] VERIFICATION (hand-derived gradients vs finite differences; beam prior vs FE)\n');
Rv = val12_ai_components();
for r = Rv
    P_('  %-58s value %10.4g  ref %10.4g  err %.1e  %s\n', r.quantity, r.value, r.ref, r.err, ifelse_str(r.pass, 'PASS', 'FAIL'));
end
gc1 = Rv(1).err; gc2 = Rv(2).err;

%% 3. PI-Net: final model ------------------------------------------------------
P_('\n[3] PI-NET - physics-informed multi-task damage network\n');
D = prep_pinet_data(F, db, S, 1.0, 1, 5, 'fbg');
tr = D.tr; va = D.va; te = D.te; htr = tr & y == 1;
P_('  inputs: %d (117 healthy-whitened features + %d family novelty indices); records train/val/test %d/%d/%d\n', ...
   size(D.X,2), size(D.X,2) - 117, sum(tr), sum(va), sum(te));
wts = struct('cls',1,'face',0.3,'map',0.3,'loc',1,'ext',0.3,'rf',0.3,'phys',LAMBDA,'l2',1e-2);
opt = struct('epochs', EPOCHS, 'noise', 0.6, 'seed', 1);
[th0, arch] = pinet_init(size(D.X,2), HIDDEN, nc, mesh.ne, 1);
tic; [th, hist] = train_pinet(th0, arch, D.X(tr,:), sub(D.T, tr), D.X(va,:), sub(D.T, va), D.phys, wts, opt);
P_('  architecture %d-%s-heads(7+1+%d+1+1), %d parameters; trained %d epochs (best %d) in %.0f s\n', ...
   size(D.X,2), mat2str(HIDDEN), mesh.ne, numel(th), numel(hist.val), hist.best_epoch, toc);
out = pinet_predict(th, arch, D.X, D.phys, D.T.E);
M = eval_pinet(out, D, F, db, te, va);
Pk = nn_unpack(th, arch.shapes);
P_('  learned evidence-prior gain a = %.2f\n', Pk{end});
fig = new_fig(800, 380);
plot(hist.train, 'b-', 'LineWidth', 1.3); hold on; plot(hist.val, 'r-', 'LineWidth', 1.3);
plot(hist.best_epoch, hist.val(hist.best_epoch), 'ko', 'MarkerSize', 9, 'LineWidth', 2);
xlabel('epoch'); ylabel('objective (supervised + \lambda physics + L2)'); grid on;
legend({'training (noisy inputs)', 'validation', 'best epoch'}, 'Location', 'northeast');
title('PI-Net training history');
save_fig(fig, 'M5_training_curves');

% Milestone 4 baselines on the same test rows
Z = double(F.X); g = F.groups;
base = {'Mahalanobis: MSE index', {'msedi'}; 'Mahalanobis: all physics', {'freq','ratio','mac','uls','strain','msedi'}};
AUCb = nan(2, nc); DRb = AUCb; FAb = nan(2,1);
for b = 1:2
    c_ = ismember(g, base{b,2});
    [~, mdl] = novelty_mahal(Z(htr, c_), Z(1, c_)); sc = novelty_mahal([], Z(:, c_), mdl);
    thr = prctile_simple(sc(va & y == 1), 99); FAb(b) = mean(sc(te & y == 1) > thr);
    for c = 2:nc
        q = te & (y == c | y == 1);
        AUCb(b,c) = roc_auc(sc(q), y(q) == c); DRb(b,c) = mean(sc(te & y == c) > thr);
    end
end
P_('\n  TEST SET - detection (AUC healthy vs class; detection rate at 1 %% false alarms)\n');
P_('  %-28s', 'method'); P_(' %9s', cn{2:end}); P_('  %6s\n', 'FA');
for b = 1:2
    P_('  %-28s', base{b,1}); P_(' %9.2f', AUCb(b,2:end)); P_('  %5.1f%%\n', 100*FAb(b));
end
P_('  %-28s', 'PI-Net (1 - p_healthy)'); P_(' %9.2f', M.auc(2:end)); P_('  %5.1f%%\n', 100*M.fa);
P_('  detection rate @1%% FA:\n  %-28s', base{1,1}); P_(' %8.0f%%', 100*DRb(1,2:end)); P_('\n');
P_('  %-28s', 'PI-Net'); P_(' %8.0f%%', 100*M.dr(2:end)); P_('\n');
P_('\n  TEST SET - 7-class diagnosis: accuracy %.3f, macro-F1 %.3f (chance 0.143)\n', M.acc, M.macroF1);
tmp = [cn; num2cell(M.f1)];
P_('  per-class F1: %s\n', sprintf('%s %.2f  ', tmp{:}));
fig = new_fig(760, 640);
Cn = M.conf ./ max(sum(M.conf, 2), 1);
imagesc(Cn); colormap(flipud(gray)); caxis([0 1]); colorbar; axis square;
set(gca, 'XTick', 1:nc, 'XTickLabel', cn, 'YTick', 1:nc, 'YTickLabel', cn);
for i = 1:nc, for j = 1:nc
    text(j, i, sprintf('%.2f', Cn(i,j)), 'HorizontalAlignment', 'center', 'Color', [0.8 0.1 0.1]*(Cn(i,j) < 0.5) + [1 1 1]*(Cn(i,j) >= 0.5));
end, end
xlabel('predicted'); ylabel('true'); title(sprintf('PI-Net test confusion (row-normalised), accuracy %.2f', M.acc));
save_fig(fig, 'M5_confusion');
fig = new_fig(1000, 400);
bar([AUCb(:,2:end); M.auc(2:end)]'); hold on; plot([0.5 nc-0.5], [0.5 0.5], 'k:');
set(gca, 'XTickLabel', cn(2:end)); ylim([0.4 1]); ylabel('test AUC'); grid on;
legend({'Mahalanobis: MSE index (M4)', 'Mahalanobis: all features (M4)', 'PI-Net'}, 'Location', 'northeast');
title('Damage detection: unsupervised physics features vs physics-informed network');
save_fig(fig, 'M5_detection_auc');

% localisation, quantification, residual strength
xc = elem_centroids(mesh);
Xc = double(F.X_comp); im = find(strcmp(g, 'msedi')); is = find(strcmp(g, 'strain'));
zz = @(A) (A - mean(A(htr,:))) ./ std(A(htr,:));
zsens = zz(Xc(:, im)) + zz(abs(Xc(:, is(1:30)))) + zz(abs(Xc(:, is(31:60))));
sig = db.P(sid,17) > 0.5 & db.P(sid,15) > 0.03;
Dd = prep_pinet_data(F, db, S, 1.0, 1, 5, 'dense');
[thd, ~] = train_pinet(th0, arch, Dd.X(tr,:), sub(Dd.T, tr), Dd.X(va,:), sub(Dd.T, va), Dd.phys, wts, opt);
outd = pinet_predict(thd, arch, Dd.X, Dd.phys, Dd.T.E); Md = eval_pinet(outd, Dd, F, db, te, va);
P_('\n  TEST SET - localisation of significant damage (D > 0.5, r > 30 mm): median error / hit < 50 mm\n');
P_('  %-14s %17s %17s %17s %17s\n', 'class', 'M4 rule (60 FBG)', 'PI-Net (60 FBG)', 'M4 rule (dense)', 'PI-Net (dense)');
loc_tab = [];
for c = 2:5
    q = find(te & y == c & sig); e1 = zeros(numel(q),1); e2 = e1;
    for j = 1:numel(q)
        e1(j) = norm(localize_damage(zsens(q(j),:)', F.sensors.fbg_pos, 2) - F.y_loc(q(j),:));
        e2(j) = norm(localize_damage(double(F.map_z(q(j),:))', xc, 3) - F.y_loc(q(j),:));
    end
    r = [median(e1) mean(e1 < 0.05) M.loc_med(2,c) M.loc_hit50(2,c) median(e2) mean(e2 < 0.05) Md.loc_med(2,c) Md.loc_hit50(2,c)];
    loc_tab(end+1,:) = [c numel(q) r]; %#ok<AGROW>
    P_('  %-14s %9.0f mm %4.0f%% %9.0f mm %4.0f%% %9.0f mm %4.0f%% %9.0f mm %4.0f%%\n', cn{c}, 1e3*r(1), 100*r(2), 1e3*r(3), 100*r(4), 1e3*r(5), 100*r(6), 1e3*r(7), 100*r(8));
end
P_('  damaged face (skin + delamination classes, test): accuracy %.2f\n', ...
   mean((out.pf(te & ismember(y,[3 4 6])) > 0.5) == F.y_face(te & ismember(y,[3 4 6]))));
P_('  damage extent (sum D / 5): R2 %.2f, Spearman rho %.2f\n', M.ext_r2, M.ext_rho);
P_('  residual strength RF_min: RMSE %.2f (R2 %.2f); clearly weakened wings (n = %d): RMSE %.2f\n', ...
   M.rf_rmse, M.rf_r2, M.n_weak, M.rf_rmse_weak);
nullp = mean(mean(D.T.dXn(te,:).^2, 2));
P_('  physics consistency: decoded-vs-measured residual %.3f (no-damage prediction: %.3f)\n', M.phys, nullp);

% example maps
big = db.P(sid,15) > 0.035 & db.P(sid,17) > 0.6;
ex = [find(te & big & y == 4 & F.y_face, 1), find(te & big & y == 3 & ~F.y_face, 1), find(te & big & y == 2, 1)];
fig = new_fig(1400, 820);
for j = 1:numel(ex)
    i = ex(j);
    subplot(3,3,3*j-2); plot_field(mesh, double(F.Dmap(i,:))', sprintf('TRUE: %s, D = %.2f', cn{y(i)}, F.y_sev(i)), 'D');
    subplot(3,3,3*j-1); plot_field(mesh, D.T.E(i,:)', 'physics evidence map E (60 FBG)', 'z');
    subplot(3,3,3*j);   plot_field(mesh, (D.phys.s_scale*out.s(i)*out.q(i,:))', ...
        sprintf('PI-Net: p(%s) = %.2f', cn{y(i)}, out.pc(i, y(i))), 'D est.');
    hold on; plot(out.loc(i,1), out.loc(i,2), 'kx', 'MarkerSize', 14, 'LineWidth', 3);
    plot(F.y_loc(i,1), F.y_loc(i,2), 'wo', 'MarkerSize', 12, 'LineWidth', 2);
end
save_fig(fig, 'M5_damage_maps');
fig = new_fig(900, 420);
rfp = out.rf*D.rf_sd + D.rf_mu;
for c = 1:nc
    q = te & y == c; plot(F.y_RF(q), rfp(q), '.', 'MarkerSize', 7); hold on;
end
lim = [min(F.y_RF(te)) max(F.y_RF(te))]; plot(lim, lim, 'k--');
xlabel('true RF_{min} (FE, ultimate load)'); ylabel('PI-Net estimate'); grid on;
legend([cn, {'ideal'}], 'Location', 'eastoutside'); title('Residual strength from sensor data (test)');
save_fig(fig, 'M5_rf_parity');

%% 4. Ablations ------------------------------------------------------------------
P_('\n[4] ABLATION - evidence-map prior (test set, lambda = %.1f)\n', LAMBDA);
Dn = prep_pinet_data(F, db, S, 1.0, 1, 5, 'none');
[thn, ~] = train_pinet(th0, arch, Dn.X(tr,:), sub(Dn.T, tr), Dn.X(va,:), sub(Dn.T, va), Dn.phys, wts, opt);
Mn = eval_pinet(pinet_predict(thn, arch, Dn.X, Dn.phys, Dn.T.E), Dn, F, db, te, va);
abl = {'no prior', Mn; 'FBG evidence prior', M; 'dense evidence prior', Md};
P_('  %-22s %6s %6s %8s   %s\n', 'variant', 'acc', 'mAUC', 'loc hit', 'median loc error (mm) debond/matrix/fibre/crush');
for k = 1:3
    Mk = abl{k,2};
    P_('  %-22s %6.3f %6.3f %7.0f%%   %s\n', abl{k,1}, Mk.acc, mean(Mk.auc(2:end)), 100*mean_finite(Mk.loc_hit50(2,2:5)), sprintf('%5.0f ', 1e3*Mk.loc_med(2,2:5)));
end

%% 5. Physics loss x label efficiency -----------------------------------------
P_('\n[5] PHYSICS LOSS x LABELLED FRACTION (semi-supervised; unlabelled wings only feel the physics loss)\n');
LF = [0.05 0.2 1]; LAM = [0 0.3 1]; SEEDS = 1:2;
if QUICK, LF = [0.05 1]; LAM = [0 0.3]; SEEDS = 1; end
LE = nan(numel(LF), numel(LAM), numel(SEEDS), 5);
for a = 1:numel(LF)
    for b = 1:numel(LAM)
        for s = SEEDS
            Dl = prep_pinet_data(F, db, S, LF(a), s, 5, 'fbg');
            wl = wts; wl.phys = LAM(b);
            [tl, al] = pinet_init(size(Dl.X,2), HIDDEN, nc, mesh.ne, s);
            tl = train_pinet(tl, al, Dl.X(tr,:), sub(Dl.T, tr), Dl.X(va,:), sub(Dl.T, va), Dl.phys, wl, setfield(opt, 'seed', s));
            ol = pinet_predict(tl, al, Dl.X, Dl.phys, Dl.T.E);
            Ml = eval_pinet(ol, Dl, F, db, te, va);
            LE(a,b,s,:) = [Ml.acc, Ml.macroF1, mean(Ml.auc(2:end)), mean_finite(Ml.loc_med(2,2:5)), Ml.phys];
        end
        m_ = squeeze(mean(LE(a,b,:,:), 3));
        P_('  labelled %5.0f%%  lambda %.1f : acc %.3f  macro-F1 %.3f  mean AUC %.3f  loc %3.0f mm  phys residual %.3f\n', ...
           100*LF(a), LAM(b), m_(1), m_(2), m_(3), 1e3*m_(4), m_(5));
    end
end
fig = new_fig(1300, 380);
ttl = {'accuracy', 'macro-F1', 'mean AUC', 'median loc. error (mm)'};
for k = 1:4
    subplot(1,4,k); hold on;
    for b = 1:numel(LAM)
        v = mean(LE(:,b,:,k), 3); if k == 4, v = 1e3*v; end
        plot(100*LF, v, 'o-', 'LineWidth', 1.5);
    end
    set(gca, 'XScale', 'log'); grid on; xlabel('labelled training wings (%)'); title(ttl{k});
end
legend(arrayfun(@(l) sprintf('\\lambda_{phys} = %.1f', l), LAM, 'UniformOutput', false), 'Location', 'best');
save_fig(fig, 'M5_label_efficiency');

%% 6. Hybrid physics-AI design surrogate --------------------------------------
P_('\n[6] HYBRID PHYSICS-AI SURROGATE (design tier: layup, thicknesses, core, material, taper, sweep, root)\n');
dz = load(fullfile(root, 'results', 'doe', 'design_database.mat'), 'P', 'spec', 'freq', 'mode_type', 'mode_order', 'tipw', 'RF', 'ok');
ctxd = doe_context(dz.spec, cfg0);
Nd = numel(dz.ok); Yd = nan(Nd, 4); Pd = Yd; Fc = nan(Nd, 6);
tic;
for i = 1:Nd
    ib = find(dz.mode_type(i,:) == 1 & dz.mode_order(i,:) == 1, 1);
    it = find(dz.mode_type(i,:) == 2 & dz.mode_order(i,:) == 1, 1);
    if ~isempty(ib), Yd(i,1) = dz.freq(i,ib); end
    if ~isempty(it), Yd(i,2) = dz.freq(i,it); end
    Yd(i,3:4) = [dz.tipw(i,2,1), dz.RF(i,1)];
    cfg = sample_to_case(dz.spec, dz.P(i,:), 'healthy', cfg0);
    pr = beam_prior(cfg, ctxd.n_ult, ctxd.x_off); Pd(i,:) = pr.y; Fc(i,:) = pr.feat;
end
onames = {'f_1B', 'f_1T', 'w_tip(1g)', 'RF_TsaiWu'};
plab = {'f_{1B} (Hz)', 'f_{1T} (Hz)', 'w_{tip} 1 g (m)', 'RF Tsai-Wu'};    % plot labels (TeX)
okd = dz.ok & all(isfinite(Yd), 2) & all(isfinite(Pd), 2);
Ylog = log(Yd); Plog = log(Pd); Xr = design_inputs(dz.P, dz.spec); Xc_ = [Xr, Fc];
P_('  reduced-order model for %d designs in %.1f s (FE: %.0f s per design)\n', Nd, toc, 9);
P_('  %-10s %18s %18s\n', 'output', 'median FE/model', 'corr(log)');
for k = 1:4
    q = okd; P_('  %-10s %18.3f %18.3f\n', onames{k}, exp(median(Ylog(q,k) - Plog(q,k))), subsref(corrcoef(Ylog(q,k), Plog(q,k)), struct('type', '()', 'subs', {{1, 2}})));
end
E_ = @(Yp, q) 100*sqrt(mean((Yp(q,:) - Ylog(q,:)).^2, 1));       % RMS log error ~ % error
modes = {'physics model only', 'pure NN (design variables)', 'NN + CLT groups', 'hybrid: physics x NN correction'};
rand('state', 7); perm = randperm(Nd); perm = perm(okd(perm)); %#ok<RAND>
ted = false(Nd,1); ted(perm(1:80)) = true; pool = perm(81:end);
NTR = [20 40 80 160 240]; SS = 1:3; if QUICK, NTR = [40 160]; SS = 1; end
LC = nan(numel(NTR), 4, numel(SS), 4);
sopt = struct('hidden', [32 32], 'lr', 3e-3, 'patience', 1500, 'epochs', 4000);
if QUICK, sopt.epochs = 400; end
for a = 1:numel(NTR)
    for s = SS
        rand('state', 100 + s); pp = pool(randperm(numel(pool))); %#ok<RAND>
        trd = false(Nd,1); trd(pp(1:NTR(a))) = true;
        vad = false(Nd,1); vad(pp(NTR(a)+1:min(numel(pp), NTR(a)+max(10, round(NTR(a)/3))))) = true;
        LC(a,1,s,:) = E_(Plog, ted);
        for m = 2:4
            Yp = surrogate_select(m, Xr, Xc_, Ylog, Plog, trd, vad, sopt, s);
            LC(a,m,s,:) = E_(Yp, ted);
        end
    end
end
P_('  test RMS error (%%, 80 unseen designs, mean of %d splits) vs number of FE training designs:\n', numel(SS));
for a = 1:numel(NTR)
    P_('   n = %3d :', NTR(a));
    for m = 1:4, P_('  [%s] %s', sprintf('%d', m), sprintf('%5.1f', squeeze(mean(LC(a,m,:,:), 3)))); end
    P_('\n');
end
P_('   [1] %s  [2] %s  [3] %s  [4] %s;  outputs: %s\n', modes{:}, strjoin(onames, ', '));
fig = new_fig(1300, 360);
for k = 1:4
    subplot(1,4,k); hold on;
    for m = 1:4, plot(NTR, mean(LC(:,m,:,k), 3), 'o-', 'LineWidth', 1.5); end
    set(gca, 'XScale', 'log'); grid on; xlabel('FE training designs'); ylabel('test RMS error (%)'); title(plab{k});
    if k == 1, legend(modes, 'Location', 'northeast'); end
end
save_fig(fig, 'M5_surrogate_learning');

% extrapolation: train on thin/moderate cores and low sweep, test outside the range
names_d = {dz.spec.vars.name};
tco = dz.P(:, strcmp(names_d, 'tcore_mm')); swp = dz.P(:, strcmp(names_d, 'sweep'));
cases = {'core thickness: train <= 20 mm, test > 24 mm', tco <= 20, tco > 24;
         'sweep: train <= 12 deg, test > 15 deg', swp <= 12, swp > 15};
P_('  EXTRAPOLATION (test designs outside the training range), RMS error %%:\n');
EX = nan(2, 4, 4);
for k = 1:2
    inr = find(okd & cases{k,2}); rand('state', 11); inr = inr(randperm(numel(inr))); %#ok<RAND>
    nv = round(0.2*numel(inr));
    vad = false(Nd,1); vad(inr(1:nv)) = true; trd = false(Nd,1); trd(inr(nv+1:end)) = true;
    tex = okd & cases{k,3};
    EX(k,1,:) = E_(Plog, tex);
    for m = 2:4, EX(k,m,:) = E_(surrogate_select(m, Xr, Xc_, Ylog, Plog, trd, vad, sopt, 1), tex); end
    P_('   %s (n_train %d, n_test %d)\n', cases{k,1}, sum(trd), sum(tex));
    for m = 1:4, P_('     %-34s %s\n', modes{m}, sprintf('%7.1f', squeeze(EX(k,m,:)))); end
end
fig = new_fig(1200, 380);
for k = 1:2
    subplot(1,2,k); bar(squeeze(EX(k,:,:))'); set(gca, 'XTickLabel', {'f1B', 'f1T', 'w tip', 'RF TW'}); grid on;
    ylabel('test RMS error (%)'); title(['Extrapolation - ' cases{k,1}]);
end
legend(modes, 'Location', 'northwest');
save_fig(fig, 'M5_surrogate_extrapolation');
% final surrogate (n = 240) and parity plot
rand('state', 101); pp = pool(randperm(numel(pool))); %#ok<RAND>
trd = false(Nd,1); trd(pp(1:240)) = true; vad = false(Nd,1); vad(pp(241:end)) = true;
[Yh, sur] = surrogate_select(4, Xr, Xc_, Ylog, Plog, trd, vad, sopt, 1);
Yn = surrogate_select(2, Xr, Xc_, Ylog, Plog, trd, vad, sopt, 1);
fig = new_fig(1300, 340);
for k = 1:4
    subplot(1,4,k);
    plot(Yd(ted,k), exp(Yn(ted,k)), 'o', 'MarkerSize', 4); hold on;
    plot(Yd(ted,k), exp(Yh(ted,k)), '.', 'MarkerSize', 12);
    lim = [min(Yd(ted,k)) max(Yd(ted,k))]; plot(lim, lim, 'k--'); grid on; axis square;
    xlabel(['FE ' plab{k}]); ylabel('surrogate'); title(plab{k});
end
legend({'pure NN', 'hybrid'}, 'Location', 'northwest');
save_fig(fig, 'M5_surrogate_parity');

%% 7. Save models ---------------------------------------------------------------
pinet = struct('theta', th, 'arch', arch, 'xinfo', D.xinfo, 'phys', D.phys, 'rf_mu', D.rf_mu, 'rf_sd', D.rf_sd, ...
    'feat_mu', D.feat_mu, 'feat_sd', D.feat_sd, 'cols', D.cols, 'emode', 'fbg', 'weights', wts, 'hist', hist, ...
    'class_names', {F.class_names}, 'test_metrics', M);
save(fullfile(outdir, 'pinet_fbg.mat'), '-struct', 'pinet', '-v7');
pinet_d = struct('theta', thd, 'arch', arch, 'emode', 'dense', 'test_metrics', Md);
save(fullfile(outdir, 'pinet_dense.mat'), '-struct', 'pinet_d', '-v7');
sur.outputs = onames; sur.inputs = 'design_inputs(P, spec) + beam_prior feat + log beam_prior y';
save(fullfile(outdir, 'surrogate_hybrid.mat'), '-struct', 'sur', '-v7');
save(fullfile(outdir, 'M5_results.mat'), 'dec_tab', 'gc1', 'gc2', 'AUCb', 'DRb', 'FAb', 'loc_tab', 'LE', 'LF', 'LAM', ...
     'LC', 'NTR', 'EX', 'modes', 'onames', 'M', 'Md', 'Mn', '-v7');
P_('\n[7] FILES: results/models/pinet_fbg.mat, pinet_dense.mat, surrogate_hybrid.mat, M5_results.mat\n');
P_('\nMilestone 5 complete (%.1f min).\n', toc(t_all)/60);
fclose(fid);
