function generate_database(tier, chunk, nchunks, outdir)
%GENERATE_DATABASE  Run the high-fidelity simulations of one DOE chunk.
%
%   generate_database(tier, chunk, nchunks, outdir)
%   Samples chunk, chunk+nchunks, chunk+2 nchunks, ... of <tier>_design.mat
%   are simulated with run_sample and written to
%   <outdir>/<tier>_chunk<chunk>of<nchunks>.mat.
%   * independent chunks -> run them in parallel (several MATLAB/Octave
%     processes, or parfor) and merge with merge_database
%   * check-pointed every 25 samples and RESUMABLE: re-running the same
%     command skips samples that are already done
%   * a failing sample is recorded (ok = false, message) but never stops
%     the run
des = load(fullfile(outdir, [tier '_design.mat']));
cfg0 = wing_baseline_config();
ctx = doe_context(des.spec, cfg0);
idx = chunk:nchunks:size(des.P, 1);
fn = fullfile(outdir, sprintf('%s_chunk%dof%d.mat', tier, chunk, nchunks));
if exist(fn, 'file')
    S = load(fn); res = S.res; done = S.done;
else
    res = cell(numel(idx), 1); done = false(numel(idx), 1);
end
t0 = tic; nnew = 0;
for j = 1:numel(idx)
    if done(j), continue; end
    i = idx(j);
    [cfg, dm] = sample_to_case(des.spec, des.P(i,:), des.classes{des.cls(i)}, cfg0);
    o = run_sample(cfg, dm, ctx);
    o.index = i;
    res{j} = o; done(j) = true; nnew = nnew + 1;
    if mod(nnew, 25) == 0 || j == numel(idx)
        save(fn, 'res', 'done', 'idx', '-v7');
        fprintf('[%s chunk %d/%d] %d/%d done, %.1f s/sample\n', tier, chunk, nchunks, ...
                sum(done), numel(idx), toc(t0)/nnew);
    end
end
save(fn, 'res', 'done', 'idx', '-v7');
end
