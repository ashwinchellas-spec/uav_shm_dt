function model = apply_bc(model, type)
%APPLY_BC  Define fixed / free DOFs.
%
%   model = apply_bc(model, type)
%   type = 'clamped_root'  wing: all 5 DOFs fixed at x = 0
%          'SSSS_crossply' simply supported (SS-1, Reddy):
%                          x = 0,a : v = w = psi_y = 0
%                          y = 0,b : u = w = psi_x = 0
%          'SSSS_hard'     simply supported, u=v=w=0 on all edges and the
%                          tangential rotation fixed (used for isotropic tests)
%          'CFFF'          clamped at x = 0 (same as clamped_root)
%          numeric vector  explicit list of fixed DOFs
%   Adds model.fixed, model.free.

mesh = model.mesh;
if isnumeric(type)
    fixed = type(:);
else
    s = mesh.sets;
    switch lower(type)
        case {'clamped_root', 'cfff'}
            fixed = dof_index(s.root, 1:5);
        case 'ssss_crossply'
            fixed = [dof_index([s.root; s.tip], [2 3 5]);
                     dof_index([s.le; s.te],   [1 3 4])];
        case 'ssss_hard'
            fixed = [dof_index([s.root; s.tip], [1 2 3 5]);
                     dof_index([s.le; s.te],   [1 2 3 4])];
        otherwise
            error('apply_bc: unknown type "%s"', type);
    end
end
fixed = unique(fixed);
free = setdiff((1:model.ndof)', fixed);
model.fixed = fixed; model.free = free;
model.bc = type;
end
