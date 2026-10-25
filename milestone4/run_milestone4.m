%RUN_MILESTONE4  Step 5: data processing & physics-informed feature engineering.
%
%   1. reference (baseline) state of the twin and the SHM sensor network
%   2. data cleaning / sanity checks of the Milestone 3 database
%   3. virtual measurements: R noisy realisations per wing (sensor noise,
%      gauge error, load uncertainty, modal-identification error)
%   4. physics-informed features + dense damage-index maps
%   5. leakage-free stratified split (by wing), EOV (temperature)
%      compensation and standardisation fitted on TRAINING data only
%   6. feature evaluation before any deep learning:
%        Fisher ratios, Mahalanobis novelty detection (AUC, detection rate
%        at 1 % false alarms) per feature set and damage class,
%        localisation error, severity correlation
%   7. ML-ready datasets for Milestone 5 (SHM tier + design tier)
%
%   >> setup_paths; run_milestone4

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
R_REAL = 3;                     % noisy realisations per simulated wing
SEED   = 42;
outdir = fullfile(root, 'results', 'features');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(fullfile(root, 'results', 'data', 'M4_summary.txt'), 'w');
P_ = @(varargin) logprintf(fid, varargin{:});
P_('=====================================================================\n');
P_(' MILESTONE 4 - Data processing & physics-informed feature engineering\n');
P_('=====================================================================\n');

db = load(fullfile(root, 'results', 'doe', 'shm_database.mat'));
N = numel(db.cls); nc = numel(db.classes);
cfg0 = wing_baseline_config();
ctx = doe_context(db.spec, cfg0);

%% 1. Reference state & sensor network -------------------------------------
ref = reference_state(db.spec, cfg0, ctx);
sens = sensor_layout(ref, 12, [0.2 0.5 0.8], 10);
mesh = ref.mesh; xc = elem_centroids(mesh);
rdb.P = ref.p; rdb.var_names = db.var_names; rdb.mode_type = ref.o.mode_type; rdb.freq = ref.o.freq;
rdb.Wmodes = ref.o.Wmodes; rdb.Emodes = ref.o.Emodes; rdb.eps_bot = ref.o.eps_bot; rdb.eps_top = ref.o.eps_top;
z0 = struct('T',0,'f',0,'phi',0,'ms',0,'eps',0,'gain',0,'load',0,'field',0);
rmeas = simulate_measurement(rdb, 1, sens, ref, z0);
P_('\n[1] REFERENCE STATE AND SENSOR NETWORK\n');
P_('  baseline: nominal properties, T = 20 C, k* = %.0f, healthy\n', ref.p(strcmp(db.var_names,'kroot')));
P_('  observed (flexural) modes: %s, f = %s Hz\n', strjoin(ref.labels, ' '), sprintf('%.2f ', rmeas.f));
P_('  FBG: %d fibres x %d gratings = %d strain sensors (bottom skin, eps_x) at y/c = %s\n', ...
   numel(sens.fbg_lines), sens.n_grat, numel(sens.fbg_elem), sprintf('%.1f ', sens.fbg_lines));
P_('  accelerometers: %d (EfI from %d grid candidates), Fisher cond = %.1f, max off-diag autoMAC = %.3f\n', ...
   numel(sens.acc_node), numel(sens.acc_cand), sens.fisher_cond, sens.automac_offdiag);
P_('  accelerometer positions [x y] (m): %s\n', sprintf('(%.2f,%.2f) ', sens.acc_pos'));
fig = new_fig(1100, 420); hold on;
patch('Faces', mesh.elems, 'Vertices', mesh.nodes, 'FaceColor', [0.93 0.95 1], 'EdgeColor', [0.8 0.8 0.9]);
for a = 1:numel(sens.fbg_lines)
    ii = (a-1)*sens.n_grat + (1:sens.n_grat);
    plot(sens.fbg_pos(ii,1), sens.fbg_pos(ii,2), 'r.-', 'MarkerSize', 16, 'LineWidth', 1);
end
plot(mesh.nodes(sens.acc_cand,1), mesh.nodes(sens.acc_cand,2), 'o', 'Color', [0.6 0.6 0.6]);
plot(sens.acc_pos(:,1), sens.acc_pos(:,2), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(0.02, 0.5*cfg0.geom.c_root, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
axis equal; set(gca, 'YDir', 'reverse'); xlim([-0.03 1.03]); ylim([-0.02 0.27]);
xlabel('span x [m]'); ylabel('chord y [m]');
title('Sensor network: FBG gratings (red), accelerometers by EfI (blue), candidates (grey), temperature (star)');
save_fig(fig, 'M4_sensor_layout');

%% 2. Data cleaning --------------------------------------------------------------
P_('\n[2] DATA CLEANING AND SANITY CHECKS\n');
bad = ~db.ok | any(~isfinite(db.freq), 2) | any(~isfinite(db.RF(:,1:5)), 2);
mono = all(diff(db.freq, 1, 2) > 0, 2);
[~, iu] = unique(round(db.P*1e9), 'rows');
P_('  failed / non-finite samples: %d;  frequencies not ascending: %d;  duplicate inputs: %d\n', ...
   sum(bad), sum(~mono), N - numel(iu));
keep = ~bad & mono;
good = find(keep);
if any(~keep)
    P_('  WARNING: %d sample(s) EXCLUDED (failed or non-finite): %s\n', sum(~keep), sprintf('%d ', find(~keep)));
    P_('           run  repair_database(''shm'', fullfile(root,''results'',''doe''))  to re-simulate them.\n');
end
h0 = db.cls == 1 & keep;
zf = abs(db.freq(h0,1:6) - mean(db.freq(h0,1:6))) ./ std(db.freq(h0,1:6));
P_('  healthy outliers |z| > 5 in f1..f6: %d;  RF range %.2f .. %.2f;  mass range %.3f .. %.3f kg\n', ...
   sum(any(zf > 5, 2)), min(db.RF(keep,5)), max(db.RF(keep,5)), min(db.mass(keep)), max(db.mass(keep)));

%% 3-4. Virtual measurements and features --------------------------------------
fx = fullfile(outdir, 'raw_features.mat');
% cache fingerprint: rebuilt automatically whenever the database changes
dbsig = [N, numel(good), sum(db.freq(good,1)), sum(db.RF(good,5)), R_REAL, SEED];
redo = ~exist(fx, 'file');
if ~redo
    vars = whos('-file', fx);
    redo = ~any(strcmp({vars.name}, 'dbsig'));
    if ~redo
        C = load(fx, 'dbsig');
        redo = numel(C.dbsig) ~= numel(dbsig) || any(abs(C.dbsig - dbsig) > 1e-9*abs(dbsig));
    end
end
if redo, P_('  (database changed or no feature cache -> features are recomputed)\n'); end
if redo
    randn('state', SEED); %#ok<RAND>
    nrow = numel(good)*R_REAL;
    [x1, names, groups] = extract_features(rmeas, rmeas);
    X = nan(nrow, numel(x1)); smap = zeros(nrow, mesh.ne, 'single'); mmap = smap;
    sid = zeros(nrow,1); rid = sid; pmac = nan(nrow, 5); Tm = nan(nrow,1);
    t0 = tic; k = 0;
    for i = good(:)'
        for rr = 1:R_REAL
            k = k + 1;
            m = simulate_measurement(db, i, sens, ref);
            [X(k,:), ~, ~, mp] = extract_features(m, rmeas);
            smap(k,:) = mp.strain'; mmap(k,:) = mp.msedi';
            sid(k) = i; rid(k) = rr; pmac(k,:) = m.pair_mac; Tm(k) = m.T;
        end
        if mod(i, 700) == 0, fprintf('  features: %d / %d wings (%.0f s)\n', i, N, toc(t0)); end
    end
    save(fx, 'X', 'names', 'groups', 'smap', 'mmap', 'sid', 'rid', 'pmac', 'Tm', 'dbsig', '-v7');
end
S = load(fx);
X = S.X; names = S.names; groups = S.groups; sid = S.sid; Tm = S.Tm;
nrow = size(X,1);
P_('\n[3-4] VIRTUAL MEASUREMENTS AND FEATURES\n');
P_('  %d wings x %d noisy realisations = %d measurement records, %d features each\n', numel(unique(sid)), R_REAL, nrow, numel(names));
P_('  noise model: T 0.5 C, frequencies 0.05 %%, mode shapes 1 %%, modal strain 2 %%,\n');
P_('               FBG 1 ue + 0.5 %% gain error, tip-test load 1 %%, distributed fibre 2 ue\n');
P_('  mode pairing by MAC: min paired MAC over all records = %.3f (all modes identified)\n', min(S.pmac(:)));
gn = unique(groups, 'stable');
for g = 1:numel(gn), P_('    group %-7s %3d features\n', gn{g}, sum(strcmp(groups, gn{g}))); end

%% 5. Split, EOV compensation, standardisation ------------------------------
split_s = split_by_sample(db.cls, [0.70 0.15 0.15], SEED);
split = split_s(sid);
ycls = db.cls(sid);
tr = split == 1; va = split == 2; te = split == 3;
htr = tr & ycls == 1;
P_('\n[5] SPLIT, EOV COMPENSATION, STANDARDISATION\n');
P_('  split by WING (all realisations of a wing in one subset): train %d / val %d / test %d wings\n', ...
   sum(split_s==1), sum(split_s==2), sum(split_s==3));
iT = strcmp(names, 'T');
[Xc, eov] = eov_compensate(X, Tm, htr);
Xc(:, iT) = X(:, iT);
mu = mean(Xc(tr,:), 1); sd = std(Xc(tr,:), 0, 1); sd(sd == 0) = 1;
Xz = (Xc - mu) ./ sd;
% temperature sensitivity removed (example: f1 and f2/f1)
i1 = find(strcmp(names, 'f1')); i2 = find(strcmp(names, 'f2/f1'));
cr = @(a) corrcoef(a, Tm(htr));
c1 = cr(X(htr,i1)); c1c = cr(Xc(htr,i1)); c2 = cr(X(htr,i2)); c2c = cr(Xc(htr,i2));
P_('  correlation with temperature (healthy train): f1 %.3f -> %.3f;  f2/f1 %.3f -> %.3f\n', ...
   c1(1,2), c1c(1,2), c2(1,2), c2c(1,2));
P_('  temperature slope of f1: %.4f Hz/C (%.3f %%/C)\n', eov.coef(2,i1), 100*eov.coef(2,i1)/mean(X(htr,i1)));
fig = new_fig(1100, 420);
subplot(1,2,1);
plot(Tm(htr), X(htr,i1), 'b.', Tm(htr), X(htr,i1) - eov.coef(2,i1)*(Tm(htr) - 20), 'r.');
grid on; xlabel('measured temperature [C]'); ylabel('f_1 [Hz]');
legend('raw', 'temperature-compensated (shifted to 20 C)', 'Location', 'southwest');
title('EOV compensation by regression (healthy training wings)');
subplot(1,2,2);
ig = find(strcmp(groups, 'msedi'));
hh = plot(1:numel(ig), Xc(htr & sid == min(sid(htr)), ig)', 'b.-'); hold on;
id = find(tr & ycls == 2 & db.P(sid, 17) > 0.7, 1);
hd = plot(1:numel(ig), Xc(id, ig), 'r.-', 'LineWidth', 1.5);
grid on; xlabel('FBG grating no. (1-10 fibre 1, 11-20 fibre 2, 21-30 fibre 3)'); ylabel('log MSE damage index');
legend([hh(1) hd], {'healthy wing (3 noisy repeats)', sprintf('debond at x = %.2f m, D = %.2f', ...
       db.P(sid(id),13)*cfg0.geom.span, db.P(sid(id), 17))}, 'Location', 'northwest');
title('Sensor-level MSE damage index');
save_fig(fig, 'M4_eov_msedi');

%% 6. Feature evaluation ----------------------------------------------------------
P_('\n[6] FEATURE EVALUATION (before any deep learning)\n');
sets = {'freq raw',             {'freq'},                                  false;
        'freq T-compensated',   {'freq'},                                  true;
        'freq + ratios (comp.)',{'freq','ratio'},                          true;
        'mode shapes (MAC+ULS)',{'mac','uls'},                             true;
        'FBG static strain',    {'strain'},                                true;
        'FBG MSE damage index', {'msedi'},                                 true;
        'all physics features', {'freq','ratio','mac','uls','strain','msedi'}, true};
ns = size(sets,1);
AUC = nan(ns, nc); DR = nan(ns, nc); FA = nan(ns,1);
for s = 1:ns
    cols = ismember(groups, sets{s,2});
    if sets{s,3}, Z = Xz(:, cols); else, Z = (X(:,cols) - mean(X(tr,cols))) ./ std(X(tr,cols)); end
    [~, mdl] = novelty_mahal(Z(htr,:), Z(1,:));
    sc = novelty_mahal([], Z, mdl);
    thr = prctile_simple(sc(va & ycls == 1), 99);       % threshold on VALIDATION healthy
    FA(s) = mean(sc(te & ycls == 1) > thr);
    for c = 2:nc
        q = te & (ycls == c | ycls == 1);
        AUC(s,c) = roc_auc(sc(q), ycls(q) == c);
        DR(s,c) = mean(sc(te & ycls == c) > thr);
    end
    if s == ns, sc_all = sc; thr_all = thr; end
end
P_('  Mahalanobis novelty detection on the TEST set (threshold = 99th pct of healthy validation)\n');
P_('  AUC (1.0 = perfect, 0.5 = chance):\n  %-24s', 'feature set'); P_('%10s', db.classes{2:end}); P_('%10s\n', 'false al.');
for s = 1:ns
    P_('  %-24s', sets{s,1}); P_('%10.3f', AUC(s,2:end)); P_('%9.1f%%\n', 100*FA(s));
end
P_('  detection rate at the 1 %% false-alarm threshold:\n  %-24s', 'feature set'); P_('%10s', db.classes{2:end}); P_('\n');
for s = 1:ns
    P_('  %-24s', sets{s,1}); P_('%9.0f%%', 100*DR(s,2:end)); P_('\n');
end
fig = new_fig(1100, 520);
imagesc(AUC(:,2:end)); colormap(flipud(hot(64))); caxis([0.5 1]); cb = colorbar; ylabel(cb, 'AUC');
set(gca, 'XTick', 1:nc-1, 'XTickLabel', strrep(db.classes(2:end), '_', ' '), 'YTick', 1:ns, 'YTickLabel', sets(:,1));
for s = 1:ns, for c = 2:nc
    text(c-1, s, sprintf('%.2f', AUC(s,c)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end, end
title('Damage detection AUC per feature set and damage class (test set, Mahalanobis novelty)');
save_fig(fig, 'M4_auc_heatmap');

% severity dependence of the detection score
P_('  rank correlation of the novelty score (all features) with severity x area, test set:\n   ');
sz = db.P(sid,17) .* db.P(sid,15).^2 .* db.P(sid,16);
for c = 2:nc
    q = te & ycls == c;
    P_(' %s %.2f,', db.classes{c}, spearman_rho(sc_all(q), sz(q)));
end
P_('\n');
fig = new_fig(1100, 450); hold on; cols7 = lines(nc);
for c = 2:nc
    q = te & ycls == c;
    semilogy(db.P(sid(q),17), sc_all(q), '.', 'Color', cols7(c,:), 'MarkerSize', 8);
end
plot([0 1], thr_all*[1 1], 'k--', 'LineWidth', 1.5);
set(gca, 'YScale', 'log'); grid on; xlabel('damage severity D'); ylabel('Mahalanobis D^2 (all features)');
legend([strrep(db.classes(2:end), '_', ' '), {'1 % false-alarm threshold'}], 'Location', 'eastoutside');
title('Detection score vs severity (test set)');
save_fig(fig, 'M4_score_vs_severity');

% localisation ------------------------------------------------------------------
% Dense maps are z-scored per element against the HEALTHY TRAINING population
% (Stubbs' normalisation): systematic EOV patterns and noisy elements are
% removed; then modal (bottom skin) and static (both skins) evidence is added.
smap = double(S.smap); mmap = double(S.mmap);
zm = (mmap - mean(mmap(htr,:))) ./ std(mmap(htr,:));
zs = (smap - mean(smap(htr,:))) ./ std(smap(htr,:));
zc = zm ./ std(zm, 0, 2) + zs ./ std(zs, 0, 2);
% sensor level: MSEDI (bottom FBGs) + |strain residual| bottom and top FBGs
ib = find(strcmp(groups, 'strain'));
sB = abs(Xc(:, ib(1:30))); sT = abs(Xc(:, ib(31:60))); sM = Xc(:, ig);
zz = @(A) (A - mean(A(htr,:))) ./ std(A(htr,:));
zsens = zz(sM) + zz(sB) + zz(sT);
P_('  damage localisation, test set, error = |estimated - true centre|:\n');
P_('  %-12s %-7s %-12s %18s %12s %18s %12s\n', 'class', 'face', 'subset', 'dense map median', 'hit<50mm', 'FBG (60) median', 'hit<50mm');
locc = [2 3 4 5];
loc_tab = [];
for c = locc
    for fsel = 0:1
        for sub = 0:1
            q = find(te & ycls == c & (db.P(sid,18) > 0.5) == fsel);
            if sub, q = q(db.P(sid(q),17) > 0.5 & db.P(sid(q),15) > 0.03); end
            e1 = zeros(numel(q),1); e2 = e1;
            for j = 1:numel(q)
                t = [db.P(sid(q(j)),13)*cfg0.geom.span, db.P(sid(q(j)),14)*cfg0.geom.c_root];
                e1(j) = norm(localize_damage(zc(q(j),:)', xc, 3) - t);
                e2(j) = norm(localize_damage(zsens(q(j),:)', sens.fbg_pos, 2) - t);
            end
            sname = 'all'; if sub, sname = 'D>0.5,r>30'; end
            fname = 'bottom'; if fsel, fname = 'top'; end
            P_('  %-12s %-7s %-12s %15.0f mm %11.0f%% %15.0f mm %11.0f%%\n', db.classes{c}, fname, sname, ...
               1e3*median(e1), 100*mean(e1 < 0.05), 1e3*median(e2), 100*mean(e2 < 0.05));
            loc_tab(end+1,:) = [c fsel sub numel(q) median(e1) mean(e1<0.05) median(e2) mean(e2<0.05)]; %#ok<AGROW>
        end
    end
end
P_('  (FBG grid spacing: 100 mm span-wise, 75 mm chord-wise -> sensor-level resolution limit)\n');
P_('  (face applies to skin damages; for core damages it only labels the impact side)\n');
fig = new_fig(1200, 820);
big = db.P(sid,15) > 0.035 & db.P(sid,17) > 0.6;      % clearly significant examples
ex = [find(te & big & ycls == 4 & db.P(sid,18) > 0.5, 1), ...
      find(te & big & ycls == 3 & db.P(sid,18) < 0.5, 1), ...
      find(te & big & ycls == 2, 1)];
for j = 1:3
    i = sid(ex(j));
    subplot(3,2,2*j-1);
    plot_field(mesh, double(db.Dmap(i,:))', sprintf('TRUE: %s (%s face), D = %.2f', strrep(db.classes{db.cls(i)}, '_', ' '), ...
               ifelse_str(db.P(i,18) > 0.5, 'top', 'bottom'), db.P(i,17)), 'D');
    subplot(3,2,2*j);
    plot_field(mesh, zc(ex(j),:)', 'ESTIMATED: combined z-scored damage-index map (noisy)', 'z');
    hold on; est = localize_damage(zc(ex(j),:)', xc, 3);
    plot(est(1), est(2), 'kx', 'MarkerSize', 14, 'LineWidth', 3);
end
save_fig(fig, 'M4_localisation_examples');

% Fisher discriminant ratio of single features (healthy vs damaged, train)
dmg = tr & ycls > 1 & ycls ~= 6;
FDR = (mean(Xz(dmg,:)) - mean(Xz(htr,:))).^2 ./ (var(Xz(dmg,:)) + var(Xz(htr,:)));
[fs, o] = sort(FDR, 'descend');
P_('  top-10 single features by Fisher ratio (healthy vs damaged, delam excluded):\n   ');
tmp = [names(o(1:10)); num2cell(fs(1:10))];
P_(' %s (%.3f)', tmp{:}); P_('\n');

%% 7. ML-ready datasets ---------------------------------------------------------
P_('\n[7] ML-READY DATASETS FOR MILESTONE 5\n');
F.X_raw = X; F.X_comp = Xc; F.X = single(Xz);
F.names = names; F.groups = groups; F.sample_id = sid; F.realisation = S.rid; F.split = split;
F.split_names = {'train', 'val', 'test'};
F.y_class = ycls; F.class_names = db.classes; F.y_detect = ycls > 1;
F.y_loc = [db.P(sid,13)*cfg0.geom.span, db.P(sid,14)*cfg0.geom.c_root];
F.y_loc(ycls == 1 | ycls == 7, :) = nan;
F.y_size = db.P(sid,15); F.y_sev = db.P(sid,17); F.y_face = db.P(sid,18) > 0.5;
F.y_size(ycls == 1) = 0; F.y_sev(ycls == 1) = 0;
F.y_RF = db.RF(sid, 5);
F.Dmap = db.Dmap(sid, :);
F.map_strain = S.smap; F.map_msedi = S.mmap; F.map_z = single(zc);
F.map_healthy_stats = struct('mu_m', mean(mmap(htr,:)), 'sd_m', std(mmap(htr,:)), 'mu_s', mean(smap(htr,:)), 'sd_s', std(smap(htr,:)));
F.loc_table = loc_tab;
F.map_grid = [mesh.nx mesh.ny];                     % element maps reshape to (nx x ny)
F.norm_mu = mu; F.norm_sd = sd; F.eov_coef = eov.coef;
F.sensors = sens; F.reference = rmeas; F.T_meas = Tm;
F.notes = ['X = standardised, temperature-compensated features (fit on training data). ', ...
           'Rows = noisy measurement records; split is by wing. Maps are per element (40 x 10).'];
save(fullfile(outdir, 'shm_features.mat'), '-struct', 'F', '-v7');
fcsv = fopen(fullfile(outdir, 'shm_features.csv'), 'w');
fprintf(fcsv, '%s\n', strjoin([{'sample_id', 'realisation', 'split', 'class', 'x_dmg', 'y_dmg', 'severity', 'RF_min'}, names], ','));
for k = 1:nrow
    fprintf(fcsv, '%d,%d,%s,%s,%.6g,%.6g,%.6g,%.6g', sid(k), S.rid(k), F.split_names{split(k)}, ...
            db.classes{ycls(k)}, F.y_loc(k,1), F.y_loc(k,2), F.y_sev(k), F.y_RF(k));
    fprintf(fcsv, ',%.6g', Xz(k,:)); fprintf(fcsv, '\n');
end
fclose(fcsv);
P_('  shm_features.mat / .csv: %d records x %d features + labels (class, detect, location,\n', nrow, numel(names));
P_('  size, severity, face, RF) + dense maps (strain residual, MSE index: %d x %d each)\n', nrow, mesh.ne);

% design tier: inputs -> outputs paired by mode TYPE (not index)
dz = load(fullfile(root, 'results', 'doe', 'design_database.mat'));
want = [1 1; 1 2; 1 3; 2 1; 2 2; 4 1];             % [type order]: 1B 2B 3B 1T 2T 1I
wn = {'f_1B', 'f_2B', 'f_3B', 'f_1T', 'f_2T', 'f_1I'};
Yf = nan(numel(dz.cls), size(want,1));
for i = 1:numel(dz.cls)
    for j = 1:size(want,1)
        k = find(dz.mode_type(i,:) == want(j,1) & dz.mode_order(i,:) == want(j,2), 1);
        if ~isempty(k), Yf(i,j) = dz.freq(i,k); end
    end
end
Dd.X_phys = dz.P; Dd.x_names = dz.var_names;
Dd.Y = [Yf, squeeze(dz.tipw(:,2,1)), squeeze(dz.tipw(:,2,2)), dz.twist(:,1), dz.RF(:,5), dz.mass];
Dd.y_names = [wn, {'w_tip_1g', 'w_tip_50N', 'twist_1g', 'RF_min', 'mass'}];
Dd.split = split_by_sample(ones(numel(dz.cls),1), [0.70 0.15 0.15], SEED);
Dd.x_mu = mean(dz.U(Dd.split==1,:)); Dd.x_sd = std(dz.U(Dd.split==1,:));
Dd.X = (dz.U - Dd.x_mu) ./ Dd.x_sd;
Dd.Wmodes = dz.Wmodes; Dd.mode_type = dz.mode_type; Dd.mode_order = dz.mode_order; Dd.freq = dz.freq;
save(fullfile(outdir, 'design_io.mat'), '-struct', 'Dd', '-v7');
P_('  design_io.mat: %d designs, %d inputs -> %d outputs; frequencies paired by mode TYPE\n', ...
   numel(dz.cls), size(Dd.X,2), size(Dd.Y,2));
tmp = [wn; num2cell(100*mean(isfinite(Yf),1))];
P_('  availability of each typed frequency within the first 8 modes: %s\n', sprintf('%s %.0f%%  ', tmp{:}));
P_('\nMilestone 4 complete.\n');
fclose(fid);
