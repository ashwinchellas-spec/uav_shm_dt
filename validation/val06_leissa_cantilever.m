function R = val06_leissa_cantilever()
%VAL06_LEISSA_CANTILEVER  Square isotropic cantilever plate (CFFF), nu = 0.3.
%   lambda = omega a^2 sqrt(rho h / D), thin plate (a/h = 1000).
%   Reference: Leissa, Vibration of Plates, NASA SP-160 (1969):
%     3.4917, 8.5246, 21.429, 27.331, 31.111
%   This is the closest classical analogue of the clamped wing.
lit = [3.4917 8.5246 21.429 27.331 31.111];
m = material_library('iso_unit'); a = 1; h = 1e-3;
L.mat = m; L.theta = 0; L.t = h;
sec = laminate_section(L, 5/6);
D = m.E1*h^3/(12*(1 - m.nu12^2));
ns = [8 16 32];
lam = zeros(numel(ns), 5);
for i = 1:numel(ns)
    mesh = rect_mesh(a, a, ns(i), ns(i));
    model = apply_bc(assemble_system(mesh, sec, struct('shear','hr')), 'clamped_root');
    [~, ~, om] = solve_modal(model, 5);
    lam(i,:) = om' * a^2 * sqrt(m.rho*h/D);
end
R = [];
for k = 1:5
    R = [R, vres('V06 Leissa CFFF', sprintf('lambda_%d (32x32)', k), lam(end,k), lit(k), 1e-2)]; %#ok<AGROW>
end
fig = new_fig(700, 450);
semilogy(ns, abs(lam./lit - 1)*100, 'o-', 'LineWidth', 1.5); grid on;
xlabel('elements per side'); ylabel('|error| vs Leissa [%]');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5', 'Location', 'northeast');
title('V06 cantilever square plate: convergence to Leissa (1969)');
save_fig(fig, 'V06_leissa');
end
