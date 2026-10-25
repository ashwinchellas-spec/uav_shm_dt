function R = val02_patch_test()
%VAL02_PATCH_TEST  MacNeal-Harder distorted 5-element patch (0.24 x 0.12).
%   Boundary nodes get the exact displacements of a constant-strain state;
%   interior nodes must reproduce it exactly (consistency / convergence).
%   Membrane:  u = 1e-3 (x + y/2),  v = 1e-3 (y + x/2)
%   Bending :  w = 1e-3 (1 + x + y + x^2/2 + xy/2 + y^2/2)/2,
%              psi = -grad w   (constant curvature, zero shear)
nodes = [0 0; .24 0; .24 .12; 0 .12; .04 .02; .18 .03; .16 .08; .08 .08];
mesh.nodes = nodes; mesh.elems = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8];
mesh.nn = 8; mesh.ne = 5; mesh.sec_id = ones(5,1);
m = material_library('iso_unit');
R = [];
for h = [1e-3 2e-2]
    L.mat = m; L.theta = 0; L.t = h;
    sec = laminate_section(L, 5/6);
    model = assemble_system(mesh, sec, struct('shear', 'hr'));
    x = nodes(:,1); y = nodes(:,2);
    ex = zeros(40,1);
    ex(1:5:end) = 1e-3*(x + y/2);
    ex(2:5:end) = 1e-3*(y + x/2);
    ex(3:5:end) = 1e-3*(1 + x + y + x.^2/2 + x.*y/2 + y.^2/2)/2;
    ex(4:5:end) = -1e-3*(1 + x + y/2)/2;
    ex(5:5:end) = -1e-3*(1 + y + x/2)/2;
    model = apply_bc(model, dof_index(1:4, 1:5));
    f = model.free; c = model.fixed;
    d = ex; d(f) = -model.K(f,f) \ (model.K(f,c) * ex(c));
    err = max(abs(d(f) - ex(f))) / max(abs(ex(f)));
    R = [R, vres('V02 patch test', sprintf('HR, a/h=%.0f, max interior error', 0.24/h), ...
         err, 0, 1e-8)]; %#ok<AGROW>
    R(end).err = err; R(end).pass = err < 1e-8;
end
end
