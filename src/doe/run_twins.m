function tw = run_twins(db, n_per_class, outfile)
%RUN_TWINS  Paired "healthy twin" simulations -> pure damage effect.
%   tw = run_twins(db, n_per_class, outfile)
%   For the first n_per_class damaged samples of every class, the SAME
%   material / temperature / root values are simulated WITHOUT damage.
%   twin difference = effect of the damage alone (variability cancels).
%   tw.idx (database rows), tw.f_dmg, tw.f_h (frequencies), tw.eps_d,
%   tw.eps_h (bottom-skin eps_x under the tip load, ne x n).
cfg0 = wing_baseline_config();
ctx = doe_context(db.spec, cfg0);
tw.idx = []; tw.f_dmg = []; tw.f_h = []; tw.eps_d = []; tw.eps_h = [];
for c = 2:numel(db.classes)
    ii = find(db.cls == c & db.ok & all(isfinite(db.freq), 2)); ii = ii(1:min(n_per_class, numel(ii)));
    for i = ii'
        [cfg, ~] = sample_to_case(db.spec, db.P(i,:), 'healthy', cfg0);
        o = run_sample(cfg, [], ctx);
        tw.idx(end+1,1) = i; tw.f_dmg(end+1,:) = db.freq(i,:); tw.f_h(end+1,:) = o.freq;
        tw.eps_d(:,end+1) = double(db.eps_bot(:,1,2,i)); tw.eps_h(:,end+1) = double(o.eps_bot(:,1,2));
    end
end
if nargin > 2, save(outfile, '-struct', 'tw', '-v7'); end
end
