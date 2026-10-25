function nfix = repair_database(tier, outdir)
%REPAIR_DATABASE  Re-run only the failed / non-finite samples of a database.
%   nfix = repair_database('shm', 'results/doe')
%   Uses the same design file and the (now robust) run_sample, writes the
%   new results into <tier>_database.mat and refreshes <tier>_scalars.csv
%   via the same packing code as merge_database.  Other samples are not
%   touched, so a 2800-sample database is repaired in seconds.
bad = check_database(tier, outdir);
nfix = 0;
if isempty(bad), fprintf('nothing to repair\n'); return; end
fn = fullfile(outdir, [tier '_database.mat']);
db = load(fn);
cfg0 = wing_baseline_config();
ctx = doe_context(db.spec, cfg0);
for i = bad(:)'
    [cfg, dm] = sample_to_case(db.spec, db.P(i,:), db.classes{db.cls(i)}, cfg0);
    o = run_sample(cfg, dm, ctx);
    fprintf('  sample %4d re-run: ok = %d  %s\n', i, o.ok, o.msg);
    if ~o.ok, continue; end
    db.ok(i) = true; db.time(i) = o.time; db.msg{i} = '';
    db.freq(i,:) = o.freq; db.zeta(i,:) = o.zeta; db.mode_type(i,:) = o.mode_type;
    db.mode_order(i,:) = o.mode_order; db.mass(i) = o.mass;
    db.Wmodes(:,:,i) = o.Wmodes; db.Emodes(:,:,:,i) = o.Emodes;
    db.gam_L1(i,:) = o.gam_L1; db.gam_L2(i,:) = o.gam_L2;
    db.w_static(:,:,i) = o.w_static; db.eps_bot(:,:,:,i) = o.eps_bot; db.eps_top(:,:,:,i) = o.eps_top;
    db.tipw(i,:,:) = reshape(o.tipw, 1, 3, 2); db.twist(i,:) = o.twist; db.RF(i,:) = o.RF;
    db.crit_elem(i) = o.crit_elem; db.crit_hashin_mode(i) = o.crit_hashin_mode;
    db.gust(i,:) = o.gust; db.Dmap(i,:) = o.Dmap; db.root_factor(i,:) = o.root_factor;
    db.qc(i,:) = o.qc;
    nfix = nfix + 1;
end
save(fn, '-struct', 'db', '-v7');
write_scalars_csv(db, fullfile(outdir, [tier '_scalars.csv']));
fprintf('%d of %d samples repaired and saved to %s\n', nfix, numel(bad), fn);
end
