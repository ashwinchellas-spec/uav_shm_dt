function model = build_wing_model(cfg, sections, sec_id)
%BUILD_WING_MODEL  Mesh + sections + assembly + BCs for the wing in cfg.
%
%   model = build_wing_model(cfg)                 healthy wing
%   model = build_wing_model(cfg, sections, sec_id)
%       user-supplied section table and element->section map (used by the
%       damage module: damaged elements point to degraded sections).
%
%   Boundary conditions (cfg.fe.bc):
%     'clamped_root'  u = v = w = psi_x = psi_y = 0 at x = 0
%     'elastic_root'  u = v = w = 0 at x = 0, plus distributed rotational
%                     springs on psi_x (span-wise bending) and psi_y
%                     (chord-wise) per unit root length:
%                        k_b = k* D11 / L,   k_t = k* D22 / L
%                     k* = cfg.fe.kroot (non-dimensional root fixity;
%                     k* -> inf recovers the clamp).  cfg.fe.root_factor
%                     (one value per root node, default 1) scales the
%                     springs locally, e.g. bolt loosening = factor < 1.
mesh = wing_mesh(cfg.geom, cfg.fe.nx, cfg.fe.ny);
if nargin < 2 || isempty(sections)
    sections = laminate_section(wing_layup(cfg), cfg.lam.kshear);
end
if nargin >= 3 && ~isempty(sec_id), mesh.sec_id = sec_id(:); end
model = assemble_system(mesh, sections, struct('shear', cfg.fe.shear));
if strcmpi(cfg.fe.bc, 'elastic_root')
    root = mesh.sets.root;
    [~, o] = sort(mesh.nodes(root, 2)); root = root(o);
    yr = mesh.nodes(root, 2);
    trib = zeros(numel(root), 1);                   % tributary root length
    dy = diff(yr);
    trib(1:end-1) = trib(1:end-1) + dy/2; trib(2:end) = trib(2:end) + dy/2;
    fac = ones(numel(root), 1);
    if isfield(cfg.fe, 'root_factor') && ~isempty(cfg.fe.root_factor)
        fac = cfg.fe.root_factor(:);
    end
    D0 = sections(1).D; L = cfg.geom.span;
    kb = cfg.fe.kroot * D0(1,1) / L; kt = cfg.fe.kroot * D0(2,2) / L;
    i4 = dof_index(root, 4); i5 = dof_index(root, 5);
    n = model.ndof;
    model.K = model.K + sparse(i4, i4, kb*trib.*fac, n, n) + sparse(i5, i5, kt*trib.*fac, n, n);
    model = apply_bc(model, dof_index(root, 1:3));
    model.root = struct('nodes', root, 'kb', kb, 'kt', kt, 'factor', fac);
else
    model = apply_bc(model, cfg.fe.bc);
end
model.cfg = cfg;
end
