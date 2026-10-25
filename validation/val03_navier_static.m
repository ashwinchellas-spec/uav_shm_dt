function R = val03_navier_static()
%VAL03_NAVIER_STATIC  SS-1 cross-ply [0/90/90/0] square plate, sinusoidal load.
%   Material E1/E2 = 25, G12 = G13 = 0.5 E2, G23 = 0.2 E2, nu12 = 0.25, k = 5/6.
%   Non-dimensional:  w_bar  = 100 w(a/2,b/2) E2 h^3 / (q0 a^4)
%                     sx_bar = sigma_x(a/2,b/2,h/2) h^2 / (q0 a^2)
%   References: own closed-form Navier FSDT (navier_fsdt.m) and the
%   published FSDT values (Reddy 2004, Mechanics of Laminated Composite
%   Plates and Shells): w_bar = 1.7100 (a/h=4), 0.6627 (10), 0.4912 (20),
%   0.4337 (100); CLT 0.4313.  Mesh 24x24 HR (with a 8/16/32 convergence
%   study at a/h = 10 -> figure V03_convergence).
m = material_library('reddy_e25'); a = 1; q0 = 1;
lit = [4 1.7100; 10 0.6627; 20 0.4912; 100 0.4337];
R = [];
for k = 1:size(lit,1)
    ah = lit(k,1); h = a/ah;
    sec = xply_sec(m, h);
    nav = navier_fsdt(sec, a, a, 'static', 'sin', q0);
    sc  = m.E2*h^3/(q0*a^4)*100;
    wc = fe_solve(sec, a, 24, h);          % even mesh: node at the centre
    [~, sx] = fe_solve(sec, a, 25, h);     % odd mesh: element centroid at the centre
    fn = nav.fields(a/2, a/2, h/2);
    R = [R, vres('V03 Navier static', sprintf('w_bar a/h=%d FE vs Navier', ah), wc*sc, nav.wc*sc, 5e-3), ...
            vres('V03 Navier static', sprintf('w_bar a/h=%d Navier vs Reddy', ah), nav.wc*sc, lit(k,2), 1e-3), ...
            vres('V03 Navier static', sprintf('sx_bar a/h=%d FE vs Navier', ah), sx*h^2/a^2, fn.sig(1)*h^2/a^2, 2e-2)]; %#ok<AGROW>
end
% convergence study (log-log error vs element size) for a/h = 10
h = 0.1; sec = xply_sec(m, h); nav = navier_fsdt(sec, a, a, 'static', 'sin', q0);
ns = [4 8 16 32 64]; err = zeros(size(ns)); errs = err;
fn = nav.fields(a/2, a/2, h/2);
for i = 1:numel(ns)
    [wc, sx] = fe_solve(sec, a, ns(i)+1, h);     % odd -> element at centre for stress
    err(i)  = abs(wc/nav.wc - 1);
    errs(i) = abs(sx/fn.sig(1) - 1);
end
p = polyfit(log(1./(ns+1)), log(err), 1);
fig = new_fig(700, 500);
loglog(1./(ns+1), err, 'bo-', 1./(ns+1), errs, 'rs-', 'LineWidth', 1.5); grid on;
xlabel('element size h_e / a'); ylabel('relative error');
legend('centre deflection', 'centre stress \sigma_x (element centroid)', 'Location', 'southeast');
title(sprintf('V03 convergence, [0/90/90/0], a/h=10: deflection rate %.2f', p(1)));
save_fig(fig, 'V03_convergence');
R = [R, vres('V03 Navier static', 'deflection convergence rate (expect ~2)', p(1), 2, 0.15)];
end

function sec = xply_sec(m, h)
L.mat = [m m m m]; L.theta = [0 90 90 0]; L.t = h/4*ones(1,4);
sec = laminate_section(L, 5/6);
end

function [wc, sx] = fe_solve(sec, a, n, h)
mesh = rect_mesh(a, a, n, n);
model = apply_bc(assemble_system(mesh, sec, struct('shear','hr')), 'SSSS_crossply');
F = pressure_load(mesh, @(x,y) sin(pi*x/a).*sin(pi*y/a), 3);
d = solve_static(model, F);
[e, ~, ~, N] = locate_point(mesh, [a/2 a/2]);
wc = N * d(5*(mesh.elems(e,:)-1) + 3);
res = recover_stresses(model, d);
% sigma_x at top surface (z = h/2) of the element containing the centre;
% recover_stresses stores material-axis stresses: top ply is 0 deg so s11 = sx
sx = res.sig(1, 2, sec.nply, e);
end
