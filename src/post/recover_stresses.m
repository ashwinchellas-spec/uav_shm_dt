function res = recover_stresses(model, d)
%RECOVER_STRESSES  Element-centre resultants, ply stresses and failure indices.
%
%   res = recover_stresses(model, d)
%   For every element (at its centre):
%     res.e0, res.kap   (3 x ne) mid-plane strains / curvatures
%     res.N, res.M      (3 x ne) force / moment resultants  [N/m, N]
%     res.Q             (2 x ne) transverse shear forces    [N/m]
%     res.eps_top, res.eps_bot (3 x ne) surface strains (strain-gauge /
%                       FBG readings in later milestones)
%     res.sig           (3 x 2 x nply x ne) material-axis in-plane stresses
%                       [s11 s22 t12] at ply bottom (1) and top (2)
%     res.tau_if        (2 x (nply-1) x ne) interlaminar [t_xz t_yz] at
%                       ply interfaces (equilibrium-based profile)
%     res.tsaiwu, res.hashin, res.maxstress  (nply x ne) worst over ply
%     res.hashin_mode   (nply x ne) governing Hashin mode (1 FT,2 FC,3 MT,4 MC)
%     res.ilss_index    ((nply-1) x ne) interlaminar shear index
%     res.core_shear    (1 x ne) core shear index  tau_max / S_core (0 if no core)
%
%   In-plane ply stresses come from  sigma = Qbar (e0 + z kappa); the
%   transverse shear stresses use  tau(z) = Q f(z)  where f(z) is the
%   layer-wise equilibrium profile of shear_correction.m (FSDT alone would
%   give a constant, physically wrong tau through the thickness).

mesh = model.mesh; secs = model.sections; ne = mesh.ne;
np = max(arrayfun(@(s) s.nply, secs));
res.e0 = zeros(3,ne); res.kap = res.e0; res.N = res.e0; res.M = res.e0;
res.Q = zeros(2,ne); res.eps_top = res.e0; res.eps_bot = res.e0;
res.sig = zeros(3,2,np,ne);
res.tau_if = zeros(2,max(np-1,1),ne);
res.tsaiwu = zeros(np,ne); res.hashin = res.tsaiwu; res.maxstress = res.tsaiwu;
res.hashin_mode = res.tsaiwu;
res.ilss_index = zeros(max(np-1,1),ne);
res.core_shear = zeros(1,ne);

for e = 1:ne
    sec = secs(mesh.sec_id(e));
    de  = d(model.edofs(e,:));
    e0  = model.rec.Bm(:,:,e) * de;
    kap = model.rec.Bb(:,:,e) * de;
    Q   = model.rec.Rq(:,:,e) * de;
    res.e0(:,e) = e0; res.kap(:,e) = kap; res.Q(:,e) = Q;
    res.N(:,e) = sec.A*e0 + sec.B*kap;
    res.M(:,e) = sec.B*e0 + sec.D*kap;
    res.eps_top(:,e) = e0 + sec.h/2*kap;
    res.eps_bot(:,e) = e0 - sec.h/2*kap;
    tz = sec.tprof;
    for k = 1:sec.nply
        mat = sec.layup.mat(k);
        zz = sec.z(k:k+1);
        sm = zeros(3,2);
        for q = 1:2
            em = sec.Te{k} * (e0 + zz(q)*kap);
            sm(:,q) = sec.Qm{k} * em;
        end
        res.sig(:,:,k,e) = sm;
        if strcmp(mat.type, 'core')
            seg = 5*k-4:5*k;
            tc = max(sqrt((Q(1)*tz.fx(seg)).^2 + (Q(2)*tz.fy(seg)).^2));
            res.core_shear(e) = max(res.core_shear(e), tc / mat.S12);
            res.maxstress(k,e) = max(abs(sm(:))) / min(mat.Xc, mat.S12);
        else
            fi = failure_indices(sm, mat);
            res.tsaiwu(k,e) = max(fi.tsaiwu);
            res.maxstress(k,e) = max(fi.maxstress);
            [res.hashin(k,e), j] = max(fi.hashin_max);
            res.hashin_mode(k,e) = fi.mode(j);
        end
        if k < sec.nply
            j = 5*k;                           % top point of ply k = interface
            tif = [Q(1)*tz.fx(j); Q(2)*tz.fy(j)];
            res.tau_if(:,k,e) = tif;
            % interlaminar shear index, weaker adjacent ply (material axes ignored:
            % S23 of the weaker ply is used as interlaminar strength)
            Sil = min(sec.layup.mat(k).S23, sec.layup.mat(k+1).S23);
            res.ilss_index(k,e) = norm(tif) / Sil;
        end
    end
end
end
