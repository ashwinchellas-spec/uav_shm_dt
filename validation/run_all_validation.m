%RUN_ALL_VALIDATION  Milestone 2 verification & validation of the FSDT solver.
%
%   Runs every benchmark (val01 ... val09), prints a PASS/FAIL table and
%   writes results/data/validation_summary.txt and validation_table.md.
%   Tolerance 'inf' rows are informational (e.g. the locking element).
%
%   >> run_all_validation

clear; close all; clc;
addpath(fileparts(fileparts(mfilename('fullpath'))));   % project root, whatever the current folder
root = setup_paths();
tests = {@val01_element_rank, @val02_patch_test, @val03_navier_static, ...
         @val04_navier_modal, @val05_shear_locking, @val06_leissa_cantilever, ...
         @val07_timoshenko_beam, @val08_sandwich_shear, @val09_newmark, @val10_wing_modes_beam, @val11_modal_response};
R = [];
for i = 1:numel(tests)
    t0 = tic;
    fprintf('running %-28s ... ', func2str(tests{i}));
    r = tests{i}();
    R = [R, r]; %#ok<AGROW>
    fprintf('%d checks, %.1f s\n', numel(r), toc(t0));
end

fid = fopen(fullfile(root, 'results', 'data', 'validation_summary.txt'), 'w');
fmd = fopen(fullfile(root, 'results', 'data', 'validation_table.md'), 'w');
fprintf(fmd, '| Test | Quantity | FE / computed | Reference | Rel. error | Tol. | Result |\n');
fprintf(fmd, '|---|---|---|---|---|---|---|\n');
logprintf(fid, '\n%-22s %-55s %13s %13s %10s %8s %s\n', 'TEST', 'QUANTITY', 'VALUE', 'REFERENCE', 'REL.ERR', 'TOL', 'RESULT');
npass = 0; nreq = 0;
for k = 1:numel(R)
    r = R(k);
    if isinf(r.tol), st = 'info'; tl = '-'; else
        nreq = nreq + 1; tl = sprintf('%.1e', r.tol);
        if r.pass, st = 'PASS'; npass = npass + 1; else, st = 'FAIL'; end
    end
    logprintf(fid, '%-22s %-55s %13.6g %13.6g %10.2e %8s %s\n', r.test, r.quantity, r.value, r.ref, r.err, tl, st);
    fprintf(fmd, '| %s | %s | %.6g | %.6g | %.2e | %s | %s |\n', r.test, r.quantity, r.value, r.ref, r.err, tl, st);
end
logprintf(fid, '\n%d / %d required checks passed.\n', npass, nreq);
fclose(fid); fclose(fmd);
save(fullfile(root, 'results', 'data', 'validation_results.mat'), 'R', '-v7');
