function [lc, env] = wing_load_cases(cfg, model)
%WING_LOAD_CASES  Static design load cases for the semi-wing.
%
%   [lc, env] = wing_load_cases(cfg, model)
%   lc(k).name, .n (load factor), .F (ndof x 1), .desc, .Lsemi (lift [N])
%
%   Net load at load factor n (accelerating, root-fixed frame):
%       F = F_lift(n W / 2)  -  n g M r_z
%   The second term is inertia relief from the wing's own mass (consistent
%   mass matrix times a unit vertical rigid-body vector r_z).
%   Lift uses the Schrenk span-wise distribution and the thin-airfoil
%   chord-wise distribution (centre of pressure at quarter chord).

env = flight_envelope(cfg.ac, cfg.g);
W  = cfg.ac.mass * cfg.g;
mesh = model.mesh;
rz = zeros(model.ndof, 1); rz(3:5:end) = 1;
Mrz = model.M * rz;

% unit lift (1 N on the semi-wing); the leading-edge sqrt singularity is
% integrated with the mapped rule of pressure_load (exact to 1e-5)
p1 = lift_distribution(cfg.geom, 1, 'schrenk');
F1 = pressure_load(mesh, p1, 6, true);
F1 = F1 / sum(F1(3:5:end));

mk = @(name, n, desc) struct('name', name, 'n', n, 'desc', desc, ...
        'Lsemi', n*W/2, 'F', n*W/2*F1 - n*cfg.g*Mrz);
lc = mk('LC1_cruise_1g', 1.0, 'Steady level flight at V_C');
lc(2) = mk('LC2_manoeuvre_pos', cfg.ac.n_pos, 'Symmetric pull-up, limit manoeuvre');
lc(3) = mk('LC3_manoeuvre_neg', cfg.ac.n_neg, 'Push-over, negative limit manoeuvre');
lc(4) = mk('LC4_gust_critical', env.n_limit_pos, ['Critical gust: ' env.critical]);
lc(5) = mk('LC5_ultimate', cfg.ac.SF*env.n_limit_pos, 'Ultimate = 1.5 x critical limit load');
% ground-test load: concentrated tip load at mid-chord (lab validation case)
Ft = point_load(mesh, [cfg.geom.span, cfg.geom.span*tand(cfg.geom.sweep_le) + ...
        0.5*cfg.geom.c_root*cfg.geom.taper], 3, cfg.an.tip_test_load);
lc(6) = struct('name', 'LC6_tip_test', 'n', 0, 'desc', ...
        sprintf('Ground test: %g N tip load at mid-chord', cfg.an.tip_test_load), ...
        'Lsemi', cfg.an.tip_test_load, 'F', Ft);
end
