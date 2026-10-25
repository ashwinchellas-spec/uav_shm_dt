function F = pressure_load(mesh, pfun, ngp, sing_le)
%PRESSURE_LOAD  Consistent nodal forces of a transverse pressure p(x,y).
%
%   F = pressure_load(mesh, pfun, ngp, sing_le)
%   pfun    : handle, p = pfun(x, y) [Pa], vectorised, positive along +z
%   ngp     : Gauss points per direction (default 4)
%   sing_le : true -> the pressure has an integrable 1/sqrt singularity at
%             the leading edge (thin-airfoil loading). In the leading-edge
%             elements the chord-wise coordinate is mapped as
%                 eta = -1 + 2 t^2,   d eta = 4 t dt,   t in [0,1]
%             so  p ~ (eta+1)^(-1/2) ~ 1/t  times  4t  is smooth and Gauss
%             quadrature becomes accurate again. (default false)
%   F       : 5*nn x 1 load vector (only w-DOFs non-zero),  f_a = int N_a p dA

if nargin < 3 || isempty(ngp), ngp = 4; end
if nargin < 4, sing_le = false; end
[gp, gw] = gauss_1d(ngp);
[q1, w1] = quad_rule(gp, gw, gp, gw);
% mapped rule in eta: t = (tau+1)/2, eta = -1 + 2 t^2, d eta = 2 t d tau
t = (gp + 1)/2;
[q2, w2] = quad_rule(gp, gw, -1 + 2*t.^2, gw .* 2 .* t);

isle = false(mesh.ne, 1);
if sing_le
    onle = false(mesh.nn, 1); onle(mesh.sets.le) = true;
    isle = onle(mesh.elems(:,1)) & onle(mesh.elems(:,2));
end
F = zeros(5*mesh.nn, 1);
for e = 1:mesh.ne
    if isle(e), q = q2; w = w2; else, q = q1; w = w1; end
    nd = mesh.elems(e,:);
    xy = mesh.nodes(nd, :);
    fe = zeros(1, 4);
    for k = 1:size(q, 1)
        [N, dN] = shape_q4(q(k,1), q(k,2));
        pq = N * xy;
        fe = fe + N * pfun(pq(1), pq(2)) * det(dN * xy) * w(k);
    end
    idx = 5*(nd - 1) + 3;
    F(idx) = F(idx) + fe(:);
end
end

function [q, w] = quad_rule(px, wx, py, wy)
[X, Y] = ndgrid(px, py);
W = wx(:) * wy(:)';
q = [X(:) Y(:)]; w = W(:);
end
