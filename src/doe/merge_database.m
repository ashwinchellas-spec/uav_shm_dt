function db = merge_database(tier, nchunks, outdir)
%MERGE_DATABASE  Collect all chunks of a tier into one database file.
%
%   db = merge_database(tier, nchunks, outdir)
%   Writes <outdir>/<tier>_database.mat  (MATLAB v7, readable by Python
%   scipy.io.loadmat) and <outdir>/<tier>_scalars.csv (one row per sample:
%   inputs, class, scalar outputs - convenient for quick ML / statistics).
%
%   DATABASE SCHEMA  (N samples, nn nodes, ne elements, nm modes)
%   inputs : U (N x k) unit-cube, P (N x k) physical, var_names, units
%   labels : cls (N x 1), classes, Dmap (N x ne) element damage D,
%            root_factor (N x ny+1)
%   modal  : freq (N x nm) [Hz], zeta (N x nm), mode_type (1 B,2 T,3 C,4 in-plane),
%            mode_order, Wmodes (nn x nm x N) w-mode shapes (mass-normalised),
%            Emodes (ne x 3 x nm x N) bottom-skin modal strains [ex ey gxy]
%   static : w_static (nn x 2 x N) for L1 = net 1 g flight load, L2 = 50 N tip;
%            eps_bot, eps_top (ne x 3 x 2 x N); tipw (N x 3 x 2) LE/mid/TE;
%            twist (N x 2) [deg];  gam_L1, gam_L2 (N x nm) modal forces
%   strength: RF (N x 6) = [Tsai-Wu Hashin core ILSS min(x>=x_off) min(all)]
%            at ultimate load; crit_elem, crit_hashin_mode
%   gust   : gust (N x 5) = [DAF_design, w_tip_peak/n [m], a_tip_peak/n [m/s^2],
%            eps_root_peak/n, DAF_sharp (t_g = T1)] per unit gust load factor
%   meta   : ok, msg, time, qc (N x 2) [equilibrium error, eigen residual], mass
des = load(fullfile(outdir, [tier '_design.mat']));
N = size(des.P, 1);
R = cell(N, 1);
for c = 1:nchunks
    S = load(fullfile(outdir, sprintf('%s_chunk%dof%d.mat', tier, c, nchunks)));
    for j = 1:numel(S.idx)
        if S.done(j), R{S.idx(j)} = S.res{j}; end
    end
end
miss = find(cellfun(@isempty, R));
if ~isempty(miss), error('merge_database: %d samples missing (first %d)', numel(miss), miss(1)); end
o1 = R{1};
nm = numel(o1.freq); nn = size(o1.Wmodes,1); ne = size(o1.Emodes,1);
db.U = des.U; db.P = des.P; db.var_names = des.var_names;
db.units = {des.spec.vars.unit}; db.cls = des.cls; db.classes = des.classes;
db.spec = des.spec; db.design_metrics = des.metrics;
db.ok = false(N,1); db.time = zeros(N,1); db.msg = cell(N,1);
db.freq = nan(N,nm); db.zeta = db.freq; db.mode_type = zeros(N,nm); db.mode_order = db.mode_type;
db.mass = nan(N,1); db.Wmodes = zeros(nn,nm,N,'single'); db.Emodes = zeros(ne,3,nm,N,'single');
db.gam_L1 = nan(N,nm); db.gam_L2 = db.gam_L1;
db.w_static = zeros(nn,2,N,'single'); db.eps_bot = zeros(ne,3,2,N,'single'); db.eps_top = db.eps_bot;
db.tipw = nan(N,3,2); db.twist = nan(N,2); db.RF = nan(N,6); db.crit_elem = zeros(N,1);
db.crit_hashin_mode = zeros(N,1); db.gust = nan(N,5); db.Dmap = zeros(N,ne,'single');
db.root_factor = ones(N, numel(o1.root_factor)); db.qc = nan(N,2);
for i = 1:N
    o = R{i};
    db.ok(i) = o.ok; db.time(i) = o.time; db.msg{i} = o.msg;
    db.freq(i,:) = o.freq; db.zeta(i,:) = o.zeta; db.mode_type(i,:) = o.mode_type;
    db.mode_order(i,:) = o.mode_order; db.mass(i) = o.mass;
    db.Wmodes(:,:,i) = o.Wmodes; db.Emodes(:,:,:,i) = o.Emodes;
    db.gam_L1(i,:) = o.gam_L1; db.gam_L2(i,:) = o.gam_L2;
    db.w_static(:,:,i) = o.w_static; db.eps_bot(:,:,:,i) = o.eps_bot; db.eps_top(:,:,:,i) = o.eps_top;
    db.tipw(i,:,:) = reshape(o.tipw, 1, 3, 2); db.twist(i,:) = o.twist; db.RF(i,:) = o.RF;
    db.crit_elem(i) = o.crit_elem; db.crit_hashin_mode(i) = o.crit_hashin_mode;
    db.gust(i,:) = o.gust; db.Dmap(i,:) = o.Dmap; db.root_factor(i,:) = o.root_factor;
    db.qc(i,:) = o.qc;
end
save(fullfile(outdir, [tier '_database.mat']), '-struct', 'db', '-v7');

write_scalars_csv(db, fullfile(outdir, [tier '_scalars.csv']));
end
