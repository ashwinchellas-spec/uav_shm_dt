function ref = reference_state(spec, cfg0, ctx)
%REFERENCE_STATE  Baseline ("as-designed, healthy") digital-twin state.
%
%   ref = reference_state(spec, cfg0, ctx)
%   Nominal material, thickness and alignment, T = 20 C, nominal root
%   fixity k* = 10^3 (centre of the DOE range), no damage.  Every feature
%   is measured RELATIVE to this state: in service this is the baseline
%   the twin stores at commissioning.
%   ref.o : full run_sample output;  ref.flex : indices of the flexural
%   (B/T) modes among the first n_modes;  ref.labels : their names.
names = {spec.vars.name};
p = zeros(1, spec.k);
for j = 1:spec.k
    switch names{j}
        case {'E1','E2','G12','Gc','Ec','rho_f','rho_c','tply','tcore'}, p(j) = 1;
        case 'T', p(j) = 20;
        case 'kroot', p(j) = 10^mean([spec.vars(j).lo spec.vars(j).hi]);
        otherwise, p(j) = 0;
    end
end
[cfg, dm] = sample_to_case(spec, p, 'healthy', cfg0);
ref.o = run_sample(cfg, dm, ctx);
ref.p = p;
ref.flex = find(ref.o.mode_type ~= 4);       % z-accelerometers cannot see in-plane modes
ref.flex = ref.flex(1:5);
lt = 'BTCI';
ref.labels = arrayfun(@(k) sprintf('%d%s', ref.o.mode_order(k), lt(ref.o.mode_type(k))), ...
                      ref.flex, 'UniformOutput', false);
ref.mesh = wing_mesh(cfg0.geom, cfg0.fe.nx, cfg0.fe.ny);
end
