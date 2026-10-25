function R = val07_timoshenko_beam()
%VAL07_TIMOSHENKO_BEAM  Cantilever strip (nu = 0) under a tip line load.
%   With nu = 0 the plate strip is in exact cylindrical bending, so FSDT
%   must reproduce Timoshenko beam theory:
%       delta = P L^3 / (3 E I) + P L / (k G A),   k = 5/6
%   Thick (L/h = 5) and thin (L/h = 100) cases separate the shear term.
E = 70e9; nu = 0; G = E/2; Lb = 1; b = 0.1; P = 1000;
m = struct('E1',E,'E2',E,'G12',G,'G13',G,'G23',G,'nu12',nu,'rho',2700, ...
    'Xt',1,'Xc',1,'Yt',1,'Yc',1,'S12',1,'S23',1,'type','isotropic','name','iso_nu0');
R = [];
for Lh = [5 100]
    h = Lb/Lh; Ly.mat = m; Ly.theta = 0; Ly.t = h;
    sec = laminate_section(Ly, 5/6);
    mesh = rect_mesh(Lb, b, 40, 4);
    model = apply_bc(assemble_system(mesh, sec, struct('shear','hr')), 'clamped_root');
    F = zeros(model.ndof, 1);
    tip = mesh.sets.tip; wts = ones(numel(tip),1); wts([1 end]) = 0.5; wts = wts/sum(wts);
    F(5*(tip-1)+3) = P*wts;
    d = solve_static(model, F);
    dfe = mean(d(5*(tip-1)+3));
    I = b*h^3/12; A = b*h;
    db = P*Lb^3/(3*E*I); ds = P*Lb/(5/6*G*A);
    R = [R, vres('V07 Timoshenko', sprintf('tip deflection L/h=%g (shear part %.1f%%)', Lh, 100*ds/(db+ds)), ...
         dfe, db + ds, 5e-3)]; %#ok<AGROW>
end
end
