function bad = check_database(tier, outdir)
%CHECK_DATABASE  Report failed / non-finite samples of a DOE database.
%   bad = check_database('shm', 'results/doe')
%   Prints, for every problematic sample: index, class, the stored error
%   message and which outputs are non-finite.  Returns the sample indices.
db = load(fullfile(outdir, [tier '_database.mat']), 'ok', 'msg', 'freq', 'RF', 'cls', 'classes', 'P', 'var_names');
nonfin = any(~isfinite(db.freq), 2) | any(~isfinite(db.RF(:,1:5)), 2);
bad = find(~db.ok | nonfin);
fprintf('%s database: %d samples, %d failed (ok = false), %d with non-finite outputs\n', ...
        tier, numel(db.ok), sum(~db.ok), sum(nonfin));
for i = bad(:)'
    msg = db.msg{i}; if isempty(msg), msg = '(no error message: silent non-finite result)'; end
    fprintf('  sample %4d  class %-15s ok=%d  n_NaN_freq=%d  msg: %s\n', i, db.classes{db.cls(i)}, ...
            db.ok(i), sum(~isfinite(db.freq(i,:))), msg);
    tmp = [db.var_names; num2cell(db.P(i,:))];
    fprintf('     inputs: %s\n', sprintf('%s=%.4g ', tmp{:}));
end
end
