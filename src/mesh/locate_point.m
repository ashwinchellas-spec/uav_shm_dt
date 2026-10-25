function [e, xi, eta, N] = locate_point(mesh, p)
%LOCATE_POINT  Find the element containing point p = [x y] and its natural coords.
%   [e, xi, eta, N] = locate_point(mesh, p)  - Newton inverse isoparametric map.
%   N (1x4) are the shape-function values at p (for interpolating fields).
%   Points slightly outside the planform are snapped to the closest element.
cent = zeros(mesh.ne, 2);
for a = 1:4, cent = cent + mesh.nodes(mesh.elems(:,a), :) / 4; end
[~, order] = sort(sum((cent - p).^2, 2));
best = inf;
for k = 1:min(9, mesh.ne)
    ee = order(k);
    xy = mesh.nodes(mesh.elems(ee,:), :);
    s = [0; 0];
    for it = 1:20
        [Nn, dN] = shape_q4(s(1), s(2));
        r = (Nn * xy)' - p(:);
        s = s - (dN * xy)' \ r;
        if norm(r) < 1e-13, break; end
    end
    viol = max(abs(s)) - 1;
    if viol < 1e-9
        e = ee; xi = s(1); eta = s(2); [N, ~] = shape_q4(xi, eta); return;
    elseif viol < best
        best = viol; e = ee; xi = max(-1, min(1, s(1))); eta = max(-1, min(1, s(2)));
    end
end
[N, ~] = shape_q4(xi, eta);
end
