function [sections, sec_id, dinfo] = apply_damage(cfg, mesh, damages, nlevels)
%APPLY_DAMAGE  Continuum-damage stiffness degradation of selected elements.
%
%   [sections, sec_id, dinfo] = apply_damage(cfg, mesh, damages, nlevels)
%   damages : struct array, each with
%     .type     'skin_matrix' | 'skin_fiber' | 'debond' | 'core_crush' | 'delam'
%     .center   [x y] [m]
%     .radius   [rx ry] semi-axes of the elliptical damage zone [m]
%     .severity D in [0, 1)   (property -> (1 - D) * property)
%     .profile  'uniform' (default) | 'gaussian' (D decays to D*exp(-2) at the rim)
%     .face     'top' | 'bottom' (skin / delam damage only, default 'bottom')
%   nlevels : number of quantised severity levels (default 20); damaged
%             elements share sections by level so assembly stays fast.
%
%   Damage modes (continuum damage mechanics, Matzenmiller-Lubliner-Taylor
%   style degradation of the ply / core moduli):
%     skin_matrix : matrix cracking     -> E2, G12, G23 of the face plies
%     skin_fiber  : fibre breakage      -> E1 (and E2, G12) of the face plies
%     debond      : face / core debond  -> core G13, G23 (loss of shear transfer)
%     core_crush  : impact core crush   -> all core moduli
%     delam       : inter-ply delam.    -> G13, G23 of the two plies at the
%                                          mid-face interface (ESL approximation)
%     root_loosening: boundary damage (bolt / bond loosening at the root);
%                   handled by root_loosening_factor + elastic root BC, so
%                   it leaves the sections healthy here.
%   Strengths degrade together with the stiffness they belong to
%   (E1 -> Xt, Xc;  E2 -> Yt, Yc;  G12 -> S12;  G13/G23 -> S23 (core: S12)),
%   so the recovered failure indices measure RESIDUAL strength - the label
%   used for prognostics in later milestones.
%
%   sections(1) is always the healthy laminate; sec_id maps elements.

if nargin < 4, nlevels = 20; end
layup0 = wing_layup(cfg);
sections = laminate_section(layup0, cfg.lam.kshear);
sec_id = ones(mesh.ne, 1);
dinfo.D = zeros(mesh.ne, numel(damages));

cent = zeros(mesh.ne, 2);
for a = 1:4, cent = cent + mesh.nodes(mesh.elems(:,a), :) / 4; end

nply = numel(layup0.theta); icore = find(arrayfun(@(m) strcmp(m.type, 'core'), layup0.mat));
nf = numel(cfg.lam.face_angles);
key = zeros(mesh.ne, numel(damages));     % quantised level per damage
for j = 1:numel(damages)
    dm = damages(j);
    if ~isfield(dm, 'profile') || isempty(dm.profile), dm.profile = 'uniform'; end
    r2 = ((cent(:,1) - dm.center(1))/dm.radius(1)).^2 + ((cent(:,2) - dm.center(2))/dm.radius(2)).^2;
    inside = r2 <= 1;
    if strcmpi(dm.type, 'root_loosening'), inside(:) = false; end
    D = zeros(mesh.ne, 1);
    if strcmpi(dm.profile, 'gaussian')
        D(inside) = dm.severity * exp(-2*r2(inside));
    else
        D(inside) = dm.severity;
    end
    dinfo.D(:,j) = D;
    key(:,j) = round(D / dm.severity * nlevels) .* (D > 0);
end

[ukeys, ~, ic] = unique(key, 'rows');
for u = 1:size(ukeys, 1)
    if all(ukeys(u,:) == 0), sid = 1;
    else
        layup = layup0;
        for j = 1:numel(damages)
            if ukeys(u,j) == 0 || strcmpi(damages(j).type, 'root_loosening'), continue; end
            dm = damages(j); Dl = dm.severity * ukeys(u,j) / nlevels;
            f = 1 - Dl;
            if ~isfield(dm, 'face') || isempty(dm.face), dm.face = 'bottom'; end
            if strcmpi(dm.face, 'top'), fp = nply-nf+1:nply; else, fp = 1:nf; end
            switch lower(dm.type)
                case 'skin_matrix'
                    layup = scale(layup, fp, {'E2','G12','G23'}, f);
                case 'skin_fiber'
                    layup = scale(layup, fp, {'E1','E2','G12'}, f);
                case 'debond'
                    layup = scale(layup, icore, {'G13','G23'}, f);
                case 'core_crush'
                    layup = scale(layup, icore, {'E1','E2','G12','G13','G23'}, f);
                case 'delam'
                    mid = fp(floor(nf/2)) + [0 1];
                    layup = scale(layup, mid, {'G13','G23'}, f);
                otherwise
                    error('apply_damage: unknown type %s', dm.type);
            end
        end
        sections(end+1) = laminate_section(layup, cfg.lam.kshear); %#ok<AGROW>
        sid = numel(sections);
    end
    sec_id(ic == u) = sid;
end
dinfo.Dmax = max(dinfo.D, [], 2);
dinfo.n_damaged = sum(sec_id > 1);
end

function layup = scale(layup, plies, props, f)
for k = plies
    % strengths that follow each modulus (core shear strength is stored in S12)
    if strcmp(layup.mat(k).type, 'core')
        strg = struct('E1', {{'Xt','Xc'}}, 'E2', {{'Yt','Yc'}}, 'G12', {{'S12'}}, ...
                      'G13', {{'S12','S23'}}, 'G23', {{'S12','S23'}});
    else
        strg = struct('E1', {{'Xt','Xc'}}, 'E2', {{'Yt','Yc'}}, 'G12', {{'S12'}}, ...
                      'G13', {{'S23'}}, 'G23', {{'S23'}});
    end
    done = {};
    for p = 1:numel(props)
        layup.mat(k).(props{p}) = layup.mat(k).(props{p}) * f;
        if isfield(strg, props{p})
            for q = strg.(props{p})
                if ~any(strcmp(done, q{1}))
                    layup.mat(k).(q{1}) = layup.mat(k).(q{1}) * f; done{end+1} = q{1}; %#ok<AGROW>
                end
            end
        end
    end
end
end
