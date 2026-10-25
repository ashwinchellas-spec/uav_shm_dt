function [N, dNdxi] = shape_q4(xi, eta)
%SHAPE_Q4  Bilinear Lagrange shape functions on [-1,1]^2.
%   N (1x4), dNdxi (2x4) = [dN/dxi; dN/deta].  Node order: (-1,-1),(1,-1),(1,1),(-1,1).
xs = [-1 1 1 -1]; es = [-1 -1 1 1];
N = 0.25 * (1 + xs*xi) .* (1 + es*eta);
dNdxi = 0.25 * [xs .* (1 + es*eta); es .* (1 + xs*xi)];
end
