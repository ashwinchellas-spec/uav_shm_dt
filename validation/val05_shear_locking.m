function R = val05_shear_locking()
%VAL05_SHEAR_LOCKING  Thin-plate limit: isotropic SS plate, uniform load.
%   Centre deflection FE / Navier-FSDT for a/h = 5 ... 10^4 on a fixed
%   16x16 mesh.  A locking element tends to zero as a/h grows; a
%   locking-free element stays near 1 for every a/h.
m = material_library('iso_unit'); a = 1;
ahs = [5 10 50 100 500 1000 1e4];
forms = {'hr', 'mitc4', 'sri', 'full'};
ratio = zeros(numel(forms), numel(ahs));
mesh = rect_mesh(a, a, 16, 16);
F = pressure_load(mesh, @(x,y) ones(size(x)), 2);
[e, ~, ~, N] = locate_point(mesh, [a/2 a/2]);
for j = 1:numel(ahs)
    L.mat = m; L.theta = 0; L.t = a/ahs(j);
    sec = laminate_section(L, 5/6);
    nav = navier_fsdt(sec, a, a, 'static', 'uniform', 1, 99);
    for i = 1:numel(forms)
        model = apply_bc(assemble_system(mesh, sec, struct('shear', forms{i})), 'SSSS_crossply');
        d = solve_static(model, F);
        ratio(i,j) = N * d(5*(mesh.elems(e,:)-1) + 3) / nav.wc;
    end
end
fig = new_fig(750, 500);
semilogx(ahs, ratio', 'o-', 'LineWidth', 1.6); grid on; ylim([0 1.1]);
xlabel('slenderness a/h'); ylabel('w_{FE} / w_{Navier}');
legend('HR mixed (used)', 'MITC4', 'SRI', 'full integration', 'Location', 'southwest');
title('V05 shear-locking test, 16x16 mesh, uniform load');
save_fig(fig, 'V05_shear_locking');
R = [];
for j = [1 4 7]
    R = [R, vres('V05 shear locking', sprintf('HR   a/h=%g  w/w_Navier', ahs(j)), ratio(1,j), 1, 1e-2), ...
            vres('V05 shear locking', sprintf('full a/h=%g  w/w_Navier (locks)', ahs(j)), ratio(4,j), 1, inf)]; %#ok<AGROW>
end
end
