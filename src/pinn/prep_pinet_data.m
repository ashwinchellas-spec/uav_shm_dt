function D = prep_pinet_data(F, db, S, labfrac, seed, s_scale, emode)
%PREP_PINET_DATA  Inputs, targets and normalised physics decoder for PI-Net.
%
%   D = prep_pinet_data(F, db, S, labfrac, seed, s_scale)
%   F       : Milestone 4 feature set (shm_features.mat)
%   db      : Milestone 3 database (needs Dmap, root_factor, P)
%   S       : damage_sensitivities output
%   labfrac : fraction of TRAINING wings whose labels may be used (rest is
%             unlabelled: only the physics loss sees it).  Selection is by
%             wing and stratified by class (realisations share status).
%   s_scale : extent scale (default 5)
%   emode   : evidence-map prior 'fbg' (default), 'dense' or 'none'
%   D.X (healthy-whitened inputs, novelty_inputs.m), D.T (targets incl. the
%   evidence map T.E), D.phys, D.tr/va/te (row masks), D.rf_mu/sd
if nargin < 6 || isempty(s_scale), s_scale = 5; end
if nargin < 7, emode = 'fbg'; end
n = size(F.X,1); sid = F.sample_id; y = F.y_class;
tr = F.split == 1; va = F.split == 2; te = F.split == 3;
[X, D.xinfo] = novelty_inputs(double(F.X), F.groups, tr & y == 1);
% physics columns: f2/f1..f5/f1 and the 60 FBG strain residuals (same order as S)
cols = [find(strcmp(F.groups, 'ratio')), find(strcmp(F.groups, 'strain'))];
htr = tr & y == 1;
mu = mean(F.X_comp(htr, cols), 1); sg = std(F.X_comp(htr, cols), 0, 1);
T.dXn = (F.X_comp(:, cols) - mu) ./ sg;
phys.Jn = S.J ./ sg(:);   phys.Jrootn = S.Jroot ./ sg(:);
phys.s_scale = s_scale; phys.XY = S.elem_xy; phys.ell = 0.1;
% targets
T.y = y;
T.face = nan(n,1); sk = ismember(y, [3 4 6]); T.face(sk) = double(F.y_face(sk));
Dm = double(F.Dmap); sD = sum(Dm, 2);
T.qtrue = zeros(n, size(Dm,2)); dm = ismember(y, 2:6) & sD > 0;
T.qtrue(dm,:) = Dm(dm,:) ./ sD(dm);
T.loc = F.y_loc; T.loc(~ismember(y, 2:6), :) = nan;
T.ext = zeros(n,1); T.ext(dm) = sD(dm) / s_scale;
rt = y == 7;
T.ext(rt) = mean(1 - db.root_factor(sid(rt), :), 2) / s_scale;
D.rf_mu = mean(F.y_RF(tr)); D.rf_sd = std(F.y_RF(tr));
T.rf = (F.y_RF - D.rf_mu) / D.rf_sd;
T.E = evidence_map(F, tr & y == 1, emode, S.elem_xy);
% label mask: fraction of training wings, stratified by class
rand('state', seed); %#ok<RAND>
lab = false(n,1);
for part = {tr, va}               % same labelled fraction in training and validation
    m = part{1};
    for c = unique(y(:))'
        w = unique(sid(m & y == c)); w = w(randperm(numel(w)));
        keep = w(1:max(1, round(labfrac*numel(w))));
        lab(m & ismember(sid, keep)) = true;
    end
end
T.lab = lab;
D.X = X; D.T = T; D.phys = phys; D.tr = tr; D.va = va; D.te = te; D.cols = cols;
D.feat_mu = mu; D.feat_sd = sg; D.emode = emode;
end
