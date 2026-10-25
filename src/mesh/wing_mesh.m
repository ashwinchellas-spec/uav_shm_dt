function mesh = wing_mesh(geom, nx, ny)
%WING_MESH  Structured Q4 mesh of a (tapered, swept) wing planform.
%
%   mesh = wing_mesh(geom, nx, ny)
%   geom.span      semi-span L          [m]   (x: root -> tip)
%   geom.c_root    root chord           [m]   (y: leading -> trailing edge)
%   geom.taper     c_tip / c_root       [-]   (default 1)
%   geom.sweep_le  leading-edge sweep   [deg] (default 0)
%   nx, ny         elements along span / chord
%
%   Mapping from the unit square (xi, eta) in [0,1]^2:
%       x = xi L,  c(x) = c_root (1 - (1 - taper) xi),
%       y = x tan(sweep) + eta c(x)
%
%   mesh.nodes  (nn x 2) coordinates          mesh.elems (ne x 4) CCW
%   mesh.xi, mesh.eta  parametric coordinates of every node
%   mesh.sets.root/tip/le/te  boundary node lists
%   mesh.nx, mesh.ny, mesh.geom

if ~isfield(geom, 'taper'),    geom.taper = 1;    end
if ~isfield(geom, 'sweep_le'), geom.sweep_le = 0; end

[XI, ETA] = ndgrid(linspace(0,1,nx+1), linspace(0,1,ny+1));
xi = XI(:); eta = ETA(:);         % node id = i + (j-1)(nx+1)
x  = xi * geom.span;
c  = geom.c_root * (1 - (1 - geom.taper) * xi);
y  = x * tand(geom.sweep_le) + eta .* c;

ne = nx*ny; elems = zeros(ne, 4); e = 0;
for j = 1:ny
    for i = 1:nx
        n1 = i + (j-1)*(nx+1);
        e = e + 1;
        elems(e,:) = [n1, n1+1, n1+1+(nx+1), n1+(nx+1)];
    end
end

mesh.nodes = [x y];
mesh.elems = elems;
mesh.xi = xi; mesh.eta = eta;
mesh.nx = nx; mesh.ny = ny; mesh.geom = geom;
tol = 1e-12;
mesh.sets.root = find(abs(xi) < tol);
mesh.sets.tip  = find(abs(xi - 1) < tol);
mesh.sets.le   = find(abs(eta) < tol);
mesh.sets.te   = find(abs(eta - 1) < tol);
mesh.nn = size(mesh.nodes, 1);
mesh.ne = ne;
mesh.sec_id = ones(ne, 1);      % section (laminate) id of every element
end
