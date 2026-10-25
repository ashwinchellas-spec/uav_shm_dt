function ctx = doe_context(spec, cfg0)
%DOE_CONTEXT  Quantities shared by all samples of a DOE tier.
%   - mesh sizes, number of modes, strength read-out station
%   - Rayleigh damping (alpha, beta) fitted ONCE on the nominal wing
%     (zeta = 2 % at modes 1 and 5): damping is a material/joint property,
%     so it is not re-fitted per sample - frequency shifts then also
%     change the modal damping consistently
%   - design gust duration t_g = 2H/V at the critical gust speed
%   - ultimate load factor n_ult = SF x n_limit
%   - unit lift load vector (when every sample shares the planform)
cfg = cfg0; cfg.fe.bc = 'elastic_root';
cfg.fe.kroot = 10^mean([2.3 3.7]);
model = build_wing_model(cfg);
[~, ~, om] = solve_modal(model, 5);
[ctx.alpha, ctx.beta] = rayleigh_coeffs(om(1), om(5), spec.zeta);
env = flight_envelope(cfg0.ac, cfg0.g);
ctx.tg = 2*cfg0.an.gust_H_chords*cfg0.ac.cbar / env.critical_V;
ctx.n_ult = cfg0.ac.SF * env.n_limit_pos;
ctx.x_off = spec.x_off;
ctx.n_modes = spec.n_modes;
ctx.dt = 5e-4;
ctx.ne = model.mesh.ne; ctx.nn = model.mesh.nn; ctx.ny = model.mesh.ny;
ctx.same_geom = strcmp(spec.tier, 'shm');
if ctx.same_geom
    p1 = lift_distribution(cfg0.geom, 1, 'schrenk');
    F = pressure_load(model.mesh, p1, 4, true);
    ctx.F1unit = F / sum(F(3:5:end));
end
end
