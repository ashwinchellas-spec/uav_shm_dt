function mesh = rect_mesh(a, b, nx, ny)
%RECT_MESH  Structured Q4 mesh of the rectangle [0,a] x [0,b].
%   Thin wrapper around wing_mesh (taper 1, no sweep). Extra node sets:
%   mesh.sets.x0, xa, y0, yb.
geom.span = a; geom.c_root = b; geom.taper = 1; geom.sweep_le = 0;
mesh = wing_mesh(geom, nx, ny);
mesh.sets.x0 = mesh.sets.root; mesh.sets.xa = mesh.sets.tip;
mesh.sets.y0 = mesh.sets.le;   mesh.sets.yb = mesh.sets.te;
end
