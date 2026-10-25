function R = val08_sandwich_shear()
%VAL08_SANDWICH_SHEAR  Shear correction factor and sandwich behaviour.
%   (a) homogeneous laminate: energy method must give k = 5/6 exactly
%   (b) sandwich cantilever strip vs Allen's sandwich beam theory
%         delta = P L^3 / (3 D) + P L / S
%         D = E_f b (t_f d^2 / 2 + t_f^3 / 6),   S = b G_c d^2 / t_c
%       FSDT with energy-equivalent k must match; FSDT with k = 5/6
%       is shown for contrast (it under-predicts the core-shear deflection).
R = [];
m = material_library('T300_5208'); L.mat = m; L.theta = 0; L.t = 1e-3;
s = laminate_section(L, 'auto');
R = [R, vres('V08 shear correction', 'homogeneous k_x (energy method)', s.k_energy(1), 5/6, 1e-12)];

f = material_library('T300_5208'); c = material_library('Rohacell51WF');
f.nu12 = 0;                         % isolate cylindrical bending (no anticlastic)
tf = 0.5e-3; tc = 20e-3; d = tc + tf; Lb = 0.5; b = 0.05; P = 100;
Ly.mat = [f c f]; Ly.theta = [0 0 0]; Ly.t = [tf tc tf];
mesh = rect_mesh(Lb, b, 50, 4);
D = f.E1*b*(tf*d^2/2 + tf^3/6);
S = b*c.G13*d^2/tc;
ref = P*Lb^3/(3*D) + P*Lb/S;
for kopt = {'auto', 5/6}
    sec = laminate_section(Ly, kopt{1});
    model = apply_bc(assemble_system(mesh, sec, struct('shear','hr')), 'clamped_root');
    F = zeros(model.ndof,1); tip = mesh.sets.tip;
    wts = ones(numel(tip),1); wts([1 end]) = 0.5; F(5*(tip-1)+3) = P*wts/sum(wts);
    dd = solve_static(model, F); wt = mean(dd(5*(tip-1)+3));
    if ischar(kopt{1})
        R = [R, vres('V08 sandwich', sprintf('tip deflection, k energy (k=%.3f)', sec.kx), wt, ref, 2e-2)]; %#ok<AGROW>
    else
        R = [R, vres('V08 sandwich', 'tip deflection, k = 5/6 (wrong for sandwich)', wt, ref, inf)]; %#ok<AGROW>
    end
end
end
