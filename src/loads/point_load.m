function F = point_load(mesh, p, comp, value)
%POINT_LOAD  Consistent nodal load of a concentrated force/moment at p = [x y].
%   F = point_load(mesh, p, comp, value)
%   comp : 1=Fx 2=Fy 3=Fz 4=M_psi_x 5=M_psi_y
[e, ~, ~, N] = locate_point(mesh, p);
F = zeros(5*mesh.nn, 1);
idx = 5*(mesh.elems(e,:) - 1) + comp;
F(idx) = N(:) * value;
end
