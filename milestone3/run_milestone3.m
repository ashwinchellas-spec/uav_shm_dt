%RUN_MILESTONE3  Steps 3-4: Design of Experiments + high-fidelity simulation database.
%
%   1. sample plans: stratified maximin-LHS for the 'shm' tier (7 classes x
%      400) and the 'design' tier (400 healthy designs)       -> make_design
%   2. high-fidelity FSDT simulations of every sample        -> generate_database
%      (chunked, resumable; set NCHUNKS > 1 and USE_PARFOR = true with the
%      Parallel Computing Toolbox, or start several MATLAB sessions, each
%      calling generate_database(tier, c, NCHUNKS, outdir) for its own c)
%   3. merge into  results/doe/<tier>_database.mat  + <tier>_scalars.csv
%   4. Morris screening of manufacturing / environmental variability
%   5. paired "healthy twin" runs -> pure damage effect vs. variability
%   6. database quality checks, statistics and figures      -> M3_summary.txt
%
%   Every step is skipped if its output file already exists, so the script
%   can be re-run safely (delete results/doe/* to start from scratch).
%
%   >> setup_paths; run_milestone3

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
outdir = fullfile(root, 'results', 'doe');
NCHUNKS = 4;            % must match the chunk files if they were made elsewhere
USE_PARFOR = false;     % true -> parfor over chunks (Parallel Computing Toolbox)
N_TWIN = 60;            % damaged samples per class re-run without damage
R_MORRIS = 20;          % Morris trajectories (13 runs each)

if ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(fullfile(root, 'results', 'data', 'M3_summary.txt'), 'w');
P_ = @(varargin) logprintf(fid, varargin{:});
P_('=====================================================================\n');
P_(' MILESTONE 3 - Design of Experiments & high-fidelity simulation database\n');
P_('=====================================================================\n');

%% 1-3. Designs, simulations, merge -----------------------------------------
tiers = {'shm', 'design'};
for t = 1:2
    tier = tiers{t};
    if ~exist(fullfile(outdir, [tier '_design.mat']), 'file'), make_design(tier, outdir); end
    if ~exist(fullfile(outdir, [tier '_database.mat']), 'file')
        fprintf(['\n*** %s_database.mat not found in %s\n' ...
                 '*** If you received the pre-computed database (uav_shm_dt_M3_database.zip),\n' ...
                 '*** stop now (Ctrl+C), unzip it into that folder and run again.\n' ...
                 '*** Otherwise the %s simulations are run now (resumable, ~1-2 h for the SHM tier).\n\n'], ...
                tier, outdir, tier);
        if USE_PARFOR
            parfor c = 1:NCHUNKS, generate_database(tier, c, NCHUNKS, outdir); end
        else
            for c = 1:NCHUNKS, generate_database(tier, c, NCHUNKS, outdir); end
        end
        merge_database(tier, NCHUNKS, outdir);
    end
end
db = load(fullfile(outdir, 'shm_database.mat'));
dz = load(fullfile(outdir, 'design_database.mat'));
N = numel(db.cls); nc = numel(db.classes); nm = size(db.freq, 2);
% valid samples only (failed / non-finite ones are reported and excluded)
okr = db.ok & all(isfinite(db.freq), 2) & all(isfinite(db.RF(:,1:5)), 2);

%% A. Sample plan ------------------------------------------------------------
P_('\n[A] SAMPLE PLAN\n');
P_('  SHM tier   : %d samples = %d classes x %d, %d variables (stratified maximin-LHS)\n', ...
   N, nc, db.spec.N_per_class, db.spec.k);
P_('  design tier: %d healthy designs, %d variables\n', numel(dz.cls), dz.spec.k);
P_('  %-26s %-8s %-8s %-6s %s\n', 'variable', 'lower', 'upper', 'scale', 'group');
for j = 1:db.spec.k
    v = db.spec.vars(j);
    P_('  %-26s %-8g %-8g %-6s %s\n', [v.name ' [' v.unit ']'], v.lo, v.hi, v.scale, v.group);
end
mt = db.design_metrics;
P_('  space filling (mean over classes):  min distance LHS %.3f vs random %.3f (+%.0f %%),\n', ...
   mean(mt(:,1)), mean(mt(:,2)), 100*(mean(mt(:,1))/mean(mt(:,2)) - 1));
P_('                                      CD2 discrepancy LHS %.4f vs random %.4f (-%.0f %%)\n', ...
   mean(mt(:,3)), mean(mt(:,4)), 100*(1 - mean(mt(:,3))/mean(mt(:,4))));

cfg0 = wing_baseline_config();
mesh = wing_mesh(cfg0.geom, cfg0.fe.nx, cfg0.fe.ny);
xc = elem_centroids(mesh);
cols = lines(nc);
fig = new_fig(1200, 800);
subplot(2,2,[1 2]); hold on;
for c = 2:nc
    ii = find(db.cls == c);
    if strcmp(db.classes{c}, 'root_loosening')
        plot(zeros(size(ii)) - 0.01*c, db.P(ii,14)*cfg0.geom.c_root, '.', 'Color', cols(c,:));
    else
        plot(db.P(ii,13)*cfg0.geom.span, db.P(ii,14)*cfg0.geom.c_root, '.', 'Color', cols(c,:));
    end
end
axis equal; set(gca, 'YDir', 'reverse'); xlim([-0.08 1]); ylim([0 0.25]); grid on;
xlabel('span x [m]'); ylabel('chord y [m]');
legend(strrep(db.classes(2:end), '_', ' '), 'Location', 'eastoutside');
title('Damage centres of all damaged samples (root loosening: band centre at x < 0)');
subplot(2,2,3);
ii = find(db.cls == 2);
rand('state', 1); Ur = rand(numel(ii), 2); %#ok<RAND>
plot(Ur(:,1), Ur(:,2), 'rx', db.U(ii,13), db.U(ii,14), 'bo', 'MarkerSize', 4);
axis square; grid on; xlabel('u_{x}'); ylabel('u_{y}');
legend('random (MC)', 'maximin LHS', 'Location', 'southoutside');
title('2-D projection, debond class');
subplot(2,2,4);
hist(db.P(db.cls > 1, 17), 16); grid on;
xlabel('damage severity D'); ylabel('count'); title('Severity (all damaged samples)');
save_fig(fig, 'M3_design_coverage');

%% B. Database quality -------------------------------------------------------
P_('\n[B] DATABASE QUALITY CHECKS\n');
P_('  successful samples: %d / %d (SHM), %d / %d (design)\n', sum(db.ok), N, sum(dz.ok), numel(dz.ok));
if any(~okr)
    P_('  EXCLUDED from all statistics: %d failed / non-finite sample(s): %s\n', sum(~okr), sprintf('%d ', find(~okr)));
    P_('  -> repair with: repair_database(''shm'', fullfile(root, ''results'', ''doe''))\n');
end
P_('  run time per sample: SHM %.2f s (median), design %.2f s\n', median(db.time), median(dz.time));
P_('  max equilibrium error %.1e, max eigen-residual %.1e\n', max(db.qc(:,1)), max(db.qc(:,2)));
h = db.cls == 1 & okr;
nomseq = mode(db.mode_type(h,:)*10 + db.mode_order(h,:), 1);
seq = db.mode_type*10 + db.mode_order;
sw6 = any(seq(:,1:6) ~= nomseq(1:6), 2);
tchar = 'BTCI';
lab = @(s) sprintf('%d%s', mod(s,10), tchar(floor(s/10)));
P_('  nominal mode sequence: %s\n', strjoin(arrayfun(lab, nomseq, 'UniformOutput', false), ' '));
P_('  SHM tier: samples whose first 6 modes are re-ordered (mode switching): %.1f %% healthy, %.1f %% all\n', ...
   100*mean(sw6(h)), 100*mean(sw6(okr)));
seqd = dz.mode_type*10 + dz.mode_order;
nomd = mode(seqd, 1);
P_('  design tier: %.1f %% of designs have a first-6-mode order different from the most common one\n', ...
   100*mean(any(seqd(:,1:6) ~= nomd(1:6), 2)));
P_('  -> across designs mode order changes, so features must be paired by mode TYPE /\n');
P_('     shape (MAC tracking, Milestone 4), never by index alone\n');

%% C. Morris screening ------------------------------------------------------------
fm = fullfile(outdir, 'morris.mat');
if ~exist(fm, 'file'), run_morris(R_MORRIS, fm); end
M = load(fm);
P_('\n[C] MORRIS SCREENING (healthy wing, %d trajectories = %d runs)\n', M.r, size(M.Y,1));
P_('  mu* normalised by the largest factor of each output (1.00 = most influential):\n');
P_('  %-8s', 'factor'); P_('%8s', M.outputs{:}); P_('\n');
mun = M.S.mustar ./ max(M.S.mustar, [], 1);
for i = 1:numel(M.factors)
    P_('  %-8s', M.factors{i}); P_('%8.2f', mun(i,:)); P_('\n');
end
[~, top] = sort(M.S.mustar(:,1), 'descend');
P_('  f1 ranking: %s\n', strjoin(M.factors(top(1:5)), ' > '));
fig = new_fig(1200, 800);
outs = [1 2 6 8];
for q = 1:4
    subplot(2,2,q);
    bar(mun(:, outs(q))); hold on;
    plot(1:numel(M.factors), M.S.sigma(:,outs(q))/max(M.S.mustar(:,outs(q))), 'r.', 'MarkerSize', 14);
    set(gca, 'XTick', 1:numel(M.factors), 'XTickLabel', M.factors); grid on;
    ylabel('\mu^* (bars), \sigma (dots)  normalised');
    title(sprintf('Morris screening: %s', strrep(M.outputs{outs(q)}, '_', ' ')));
    if exist('OCTAVE_VERSION', 'builtin') == 0, xtickangle(45); end
end
save_fig(fig, 'M3_morris');

% first-order Sobol indices from the healthy LHS samples (no extra runs)
P_('  first-order Sobol indices (binned estimator, %d healthy LHS samples):\n', sum(h));
Yh = [db.freq(h,[1 2 3]), db.tipw(h,2,1)];
yn = {'f1', 'f2', 'f3', 'w_tip'};
P_('  %-8s', 'factor'); P_('%8s', yn{:}); P_('\n');
S1 = zeros(12, 4);
for q = 1:4, S1(:,q) = sobol_binned(db.U(h,1:12), Yh(:,q), 10)'; end
for i = 1:12, P_('  %-8s', db.var_names{i}); P_('%8.3f', S1(i,:)); P_('\n'); end
P_('  %-8s', 'sum'); P_('%8.3f', sum(S1, 1)); P_('\n');

%% D. Damage effect vs. environmental/manufacturing variability --------------
P_('\n[D] DAMAGE SIGNAL vs VARIABILITY\n');
ft = fullfile(outdir, 'twins.mat');
if ~exist(ft, 'file'), run_twins(db, N_TWIN, ft); end
tw = load(ft);
sig_h = std(db.freq(h,:), 0, 1) ./ mean(db.freq(h,:), 1) * 100;
P_('  healthy variability (std of f_k, %%): %s\n', sprintf('%.2f ', sig_h(1:6)));
P_('  pure damage effect from paired healthy twins (|df/f| %%, median / 90th pct, worst mode):\n');
twc = db.cls(tw.idx);
dfp = abs(tw.f_dmg./tw.f_h - 1)*100;
ratio = zeros(nc,1);
for c = 2:nc
    q = twc == c; dmax = max(dfp(q,1:6), [], 2);
    srt = sort(dmax); p90 = srt(max(1, round(0.9*numel(srt))));
    ratio(c) = median(dmax ./ max(sig_h(1:6)));
    P_('    %-15s median %.3f %%, 90th pct %.3f %%  -> median signal / EOV std = %.2f\n', ...
       db.classes{c}, median(dmax), p90, ratio(c));
end
% classical outlier test using frequencies only: any f_k outside healthy mean +/- 3 std
mu_h = mean(db.freq(h,:), 1); sd_h = std(db.freq(h,:), 0, 1);
out3 = any(abs(db.freq(:,1:6) - mu_h(1:6)) > 3*sd_h(1:6), 2);
P_('  frequency-only 3-sigma novelty test: false alarms %.1f %% (healthy), detections:', 100*mean(out3(h)));
for c = 2:nc, P_(' %s %.0f %%,', db.classes{c}, 100*mean(out3(db.cls == c & okr))); end
P_('\n');
% local strain-field signal (tip-load case, bottom-skin eps_x) vs strain variability
P_('  local strain signal (paired twins, tip load): max|d eps_x| / max|eps_x|  [%%]\n');
rel = max(abs(tw.eps_d - tw.eps_h), [], 1) ./ max(abs(tw.eps_h), [], 1) * 100;
for c = 2:nc
    q = twc == c;
    P_('    %-15s median %6.2f %%, 90th pct %6.2f %%\n', db.classes{c}, median(rel(q)), ...
       prctile_simple(rel(q), 90));
end
Eh = squeeze(double(db.eps_bot(:,1,2,h)));            % ne x n_healthy
sig = abs(mean(Eh, 2)) > 0.1*max(abs(mean(Eh, 2)));
cov_e = std(Eh(sig,:), 0, 2) ./ abs(mean(Eh(sig,:), 2)) * 100;
P_('  healthy strain variability (EOV) at the same load: median CoV %.2f %% per element\n', median(cov_e));
% EOV mostly SCALES the whole strain field (stiffness, thickness, temperature);
% remove the best-fit global scale a_n of every healthy field and look at
% what is left (the variability of the strain-field SHAPE)
Em = mean(Eh(sig,:), 2);
a = (Em' * Eh(sig,:)) / (Em' * Em);
Rs = Eh(sig,:) - Em * a;
cov_s = std(Rs, 0, 2) ./ abs(Em) * 100;
P_('  after removing one global scale factor per sample: median residual CoV %.2f %%\n', median(cov_s));
relq = zeros(nc,1);
for c = 2:nc, relq(c) = median(rel(twc == c)); end
P_('  -> raw strain changes (0.2-1.4 %% median) are below the raw EOV scatter (%.1f %%), but\n', median(cov_e));
P_('     EOV is almost a global scaling: the SHAPE of the field varies only %.2f %%, so local\n', median(cov_s));
P_('     damage patterns (debond %.1fx, fibre %.1fx, matrix %.1fx, root %.1fx that level) become\n', ...
   relq(2)/median(cov_s), relq(4)/median(cov_s), relq(3)/median(cov_s), relq(7)/median(cov_s));
P_('     visible after spatial normalisation -> key design rule for Milestone 4 features\n');
P_('  delamination inside a 0.5 mm face: no measurable global effect (|df/f| <= %.1e) -\n', ...
   max(max(abs(tw.f_dmg(twc==6,:)./tw.f_h(twc==6,:) - 1))));
P_('     a physical detectability limit of low-frequency vibration / strain SHM\n');
P_('  -> natural frequencies alone cannot detect most damage under EOV: the AI needs\n');
P_('     physics-informed features (mode-shape curvature, strain fields) - Milestones 4-5\n');
fig = new_fig(1100, 500);
subplot(1,2,1); hold on;
for c = 2:nc
    q = twc == c;
    semilogy(c + 0.25*randn(sum(q),1)*0.3, max(dfp(q,1:6), [], 2), '.', 'Color', cols(c,:), 'MarkerSize', 10);
end
plot([1.5 nc+0.5], max(sig_h(1:6))*[1 1], 'k--', 'LineWidth', 1.5);
set(gca, 'YScale', 'log', 'XTick', 2:nc, 'XTickLabel', strrep(db.classes(2:end), '_', ' ')); grid on;
ylim([1e-7 1e1]);
ylabel('pure damage effect max_k |\Delta f_k / f_k| [%]');
title('Damage signal (paired twins) vs EOV std (dashed)');
subplot(1,2,2);
ex = find(twc == 2, 1);
de = (tw.eps_d(:,ex) - tw.eps_h(:,ex))*1e6;
plot_field(mesh, de, sprintf('Debond sample %d: \\Delta\\epsilon_x (tip load), bottom skin [\\mu\\epsilon]', ...
           tw.idx(ex)), '\mu\epsilon');
save_fig(fig, 'M3_damage_vs_eov');

%% E. Output statistics per class ------------------------------------------------
P_('\n[E] OUTPUT STATISTICS PER CLASS (mean +/- std)\n');
P_('  %-15s %16s %16s %16s %14s %12s %10s\n', 'class', 'f1 [Hz]', 'f2 [Hz]', 'f3 [Hz]', ...
   'w_tip 1g [mm]', 'RF_min', 'DAF_sharp');
for c = 1:nc
    q = db.cls == c & okr;
    P_('  %-15s %8.3f+/-%5.3f %8.3f+/-%5.3f %8.3f+/-%5.3f %7.4f+/-%5.4f %6.2f+/-%4.2f %5.3f\n', db.classes{c}, ...
       mean(db.freq(q,1)), std(db.freq(q,1)), mean(db.freq(q,2)), std(db.freq(q,2)), ...
       mean(db.freq(q,3)), std(db.freq(q,3)), 1e3*mean(db.tipw(q,2,1)), 1e3*std(db.tipw(q,2,1)), ...
       mean(db.RF(q,5)), std(db.RF(q,5)), mean(db.gust(q,5)));
end
fig = new_fig(1100, 450);
subplot(1,2,1); hold on;
for c = 1:nc
    q = db.cls == c & okr;
    sev = db.P(q,17); if c == 1, sev = zeros(sum(q),1); end
    plot(sev, db.RF(q,5), '.', 'Color', cols(c,:), 'MarkerSize', 8);
end
grid on; xlabel('damage severity D'); ylabel('RF_{min} at ultimate load');
legend(strrep(db.classes, '_', ' '), 'Location', 'southwest'); title('Residual strength (prognostics label)');
subplot(1,2,2);
k3 = [1 2 3];
bm = zeros(nc, 3); bs = bm;
for c = 1:nc
    q = db.cls == c & okr;
    bm(c,:) = mean(db.freq(q,k3)./mu_h(k3) - 1, 1)*100; bs(c,:) = std(db.freq(q,k3)./mu_h(k3), 0, 1)*100;
end
bar(bm); hold on; grid on;
set(gca, 'XTick', 1:nc, 'XTickLabel', strrep(db.classes, '_', ' '));
ylabel('mean (f_k / f_k^{healthy} - 1) [%]'); legend('f1 (1B)', 'f2 (1T)', 'f3 (2B)', 'Location', 'southwest');
title('Class-mean frequency shift (bars) - small compared with EOV');
save_fig(fig, 'M3_class_statistics');

%% F. Example damage maps --------------------------------------------------------
fig = new_fig(1200, 800);
for c = 2:nc
    subplot(3,2,c-1);
    i = find(db.cls == c & okr & db.P(:,17) > 0.6, 1);
    if strcmp(db.classes{c}, 'root_loosening')
        rf = db.root_factor(i,:);
        plot(rf, linspace(0, cfg0.geom.c_root, numel(rf)), 'r.-', 'LineWidth', 1.5);
        set(gca, 'YDir', 'reverse'); grid on; xlim([0 1.05]);
        xlabel('root spring factor'); ylabel('chord y [m]');
        title(sprintf('#%d root loosening, D = %.2f', i, db.P(i,17)));
    else
        plot_field(mesh, double(db.Dmap(i,:))', sprintf('#%d %s, D = %.2f, %s face', i, ...
            strrep(db.classes{c}, '_', ' '), db.P(i,17), ifelse_str(db.P(i,18) > 0.5, 'top', 'bottom')), 'D');
    end
end
save_fig(fig, 'M3_damage_examples');

%% G. Design tier -------------------------------------------------------------
P_('\n[G] DESIGN TIER (healthy wings across the design space)\n');
im = min(3, floor(dz.P(:,8)*3) + 1);
for m = 1:3
    q = im == m;
    P_('  %-13s n = %3d   f1 = %6.2f .. %6.2f Hz   w_tip(1g) = %6.2f .. %6.2f mm   RF_min = %5.2f .. %5.2f\n', ...
       dz.spec.face_mats{m}, sum(q), min(dz.freq(q,1)), max(dz.freq(q,1)), 1e3*min(dz.tipw(q,2,1)), ...
       1e3*max(dz.tipw(q,2,1)), min(dz.RF(q,5)), max(dz.RF(q,5)));
end
S1d = zeros(dz.spec.k, 1);
S1d(:) = sobol_binned(dz.U, dz.freq(:,1), 10);
[~, o] = sort(S1d, 'descend');
P_('  f1 first-order Sobol ranking: %s\n', strjoin(arrayfun(@(i) sprintf('%s (%.2f)', ...
   dz.var_names{i}, S1d(i)), o(1:5)', 'UniformOutput', false), ', '));
fig = new_fig(1100, 450);
subplot(1,2,1); hold on; mk = 'osd';
for m = 1:3
    q = im == m; plot(dz.freq(q,1), dz.freq(q,2), mk(m), 'MarkerSize', 4);
end
grid on; xlabel('f_1 [Hz]'); ylabel('f_2 [Hz]'); legend(strrep(dz.spec.face_mats, '_', '/'), 'Location', 'northwest');
title('Design tier: first two frequencies');
subplot(1,2,2); bar(S1d(o)); set(gca, 'XTick', 1:dz.spec.k, 'XTickLabel', dz.var_names(o)); grid on;
ylabel('first-order Sobol index of f_1'); title('What drives f_1 across designs');
save_fig(fig, 'M3_design_tier');

%% H. Storage -----------------------------------------------------------------------
d1 = dir(fullfile(outdir, 'shm_database.mat')); d2 = dir(fullfile(outdir, 'design_database.mat'));
P_('\n[H] FILES\n  shm_database.mat %.0f MB, design_database.mat %.0f MB, *_scalars.csv, morris.mat, twins.mat\n', ...
   d1.bytes/1e6, d2.bytes/1e6);
P_('\nMilestone 3 complete.\n');
fclose(fid);
