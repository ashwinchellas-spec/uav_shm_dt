function des = make_design(tier, outdir)
%MAKE_DESIGN  Build and save the DOE sample plan of one tier.
%   des = make_design('shm' | 'design', outdir)
%   Stratified maximin-LHS: an independent LHS of N_per_class points is
%   drawn for every damage class, so the damage parameters (location,
%   size, severity) are space-filling WITHIN each class and every class
%   has the same size (balanced classification data).
%   Saves <outdir>/<tier>_design.mat with fields
%     U (N x k, unit cube), P (N x k, physical), cls (N x 1 class index),
%     classes, spec, metrics (maximin distance and CD2 of LHS vs random).
spec = doe_config(tier);
nc = numel(spec.classes); n = spec.N_per_class;
U = zeros(spec.N, spec.k); cls = zeros(spec.N, 1);
met = zeros(nc, 4);
for c = 1:nc
    [Uc, info] = lhs_maximin(n, spec.k, spec.seed + c, 20, 400);
    rows = (c-1)*n + (1:n);
    U(rows,:) = Uc; cls(rows) = c;
    rand('state', spec.seed + 100 + c); Ur = rand(n, spec.k); %#ok<RAND>
    met(c,:) = [info.dmin, info.dmin_random, discrepancy_cd2(Uc), discrepancy_cd2(Ur)];
end
P = doe_map(spec, U);
des.U = U; des.P = P; des.cls = cls; des.classes = spec.classes; des.spec = spec;
des.var_names = {spec.vars.name};
des.metrics = met;          % [dmin_LHS dmin_random CD2_LHS CD2_random] per class
if ~exist(outdir, 'dir'), mkdir(outdir); end
save(fullfile(outdir, [tier '_design.mat']), '-struct', 'des', '-v7');
end
