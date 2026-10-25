function layup = wing_layup(cfg)
%WING_LAYUP  Full sandwich stacking sequence (bottom -> top) from cfg.lam.
%   [face | core | mirrored face], e.g. [0/45/-45/90/C/90/-45/45/0].
%
%   Optional cfg.unc (Milestone 3: manufacturing scatter and environment)
%   multiplies the nominal properties - every field is optional:
%     .E1 .E2 .G12          face modulus scale factors          [-]
%     .Gc .Ec               core shear / normal modulus scales  [-]
%     .tply .tcore          ply / core thickness scales         [-]
%     .rho_f .rho_c         face / core density scales          [-]
%     .dtheta               global ply misalignment             [deg]
%     .T                    temperature                         [deg C]
%   Temperature (reference T0 = 20 C) softens matrix-dominated and core
%   properties linearly below Tg (typical epoxy / PMI-foam behaviour):
%     E2, G12, G13, G23, Yt, Yc, S12, S23 :  x (1 - a_m (T - T0)), a_m = 0.003 /C
%     E1, Xt, Xc (fibre dominated)        :  x (1 - a_f (T - T0)), a_f = 0.0002 /C
%     core moduli and strengths           :  x (1 - a_c (T - T0)), a_c = 0.0025 /C
fm = material_library(cfg.lam.face_mat);
cm = material_library(cfg.lam.core_mat);
fa = cfg.lam.face_angles(:)';
tp = cfg.lam.t_ply; tc = 0;
if isfield(cfg.lam, 't_core'), tc = cfg.lam.t_core; end

if isfield(cfg, 'unc') && ~isempty(cfg.unc)
    u = cfg.unc;
    g = @(f, d) getdef(u, f, d);
    fm.E1 = fm.E1*g('E1',1); fm.E2 = fm.E2*g('E2',1);
    fm.G12 = fm.G12*g('G12',1); fm.G13 = fm.G13*g('G12',1); fm.G23 = fm.G23*g('E2',1);
    cm.G12 = cm.G12*g('Gc',1); cm.G13 = cm.G13*g('Gc',1); cm.G23 = cm.G23*g('Gc',1);
    cm.E1 = cm.E1*g('Ec',1); cm.E2 = cm.E2*g('Ec',1);
    fm.rho = fm.rho*g('rho_f',1); cm.rho = cm.rho*g('rho_c',1);
    tp = tp*g('tply',1); tc = tc*g('tcore',1);
    fa = fa + g('dtheta',0);
    dT = g('T',20) - 20;
    am = 0.003; af = 0.0002; ac = 0.0025;
    sm = 1 - am*dT; sf = 1 - af*dT; sc = 1 - ac*dT;
    fm.E2 = fm.E2*sm; fm.G12 = fm.G12*sm; fm.G13 = fm.G13*sm; fm.G23 = fm.G23*sm;
    fm.Yt = fm.Yt*sm; fm.Yc = fm.Yc*sm; fm.S12 = fm.S12*sm; fm.S23 = fm.S23*sm;
    fm.E1 = fm.E1*sf; fm.Xt = fm.Xt*sf; fm.Xc = fm.Xc*sf;
    for p = {'E1','E2','G12','G13','G23','Xt','Xc','Yt','Yc','S12','S23'}
        cm.(p{1}) = cm.(p{1})*sc;
    end
end

nf = numel(fa);
mats = repmat(fm, 1, 2*nf);
t = tp * ones(1, 2*nf);
if tc > 0
    layup.mat = [mats(1:nf), cm, mats(nf+1:end)];
    layup.theta = [fa, 0, fliplr(fa)];
    layup.t = [t(1:nf), tc, t(nf+1:end)];
else
    layup.mat = mats; layup.theta = [fa, fliplr(fa)]; layup.t = t;
end
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)) && ~isnan(s.(f)), v = s.(f); else, v = d; end
end
