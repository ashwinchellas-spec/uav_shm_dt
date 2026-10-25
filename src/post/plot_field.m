function plot_field(mesh, vals, ttl, clab, deform)
%PLOT_FIELD  Contour plot of a nodal (nn x 1) or element (ne x 1) field.
%   plot_field(mesh, vals, ttl, clab)
%   plot_field(mesh, vals, ttl, clab, wz)  3-D view with nodal z = wz
%   Draws into the current axes.  LE at the top (y axis reversed).
if nargin < 4, clab = ''; end
X = mesh.nodes(:,1); Y = mesh.nodes(:,2);
if nargin >= 5 && ~isempty(deform), V = [X Y deform(:)]; else, V = [X Y]; end
vals = vals(:);
faces = mesh.elems;
if size(V,2) == 3
    % split quads into triangles (Octave's gnuplot back-end draws only
    % triangular 3-D patches; MATLAB is indifferent)
    faces = [mesh.elems(:,[1 2 3]); mesh.elems(:,[1 3 4])];
    if numel(vals) == mesh.ne, vals = [vals; vals]; end
end
if numel(vals) == mesh.nn
    patch('Faces', faces, 'Vertices', V, 'FaceVertexCData', vals, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
else
    patch('Faces', faces, 'Vertices', V, 'FaceVertexCData', vals, ...
          'FaceColor', 'flat', 'EdgeColor', 'none');
end
colormap(jet(64)); cb = colorbar; ylabel(cb, clab);
if size(V,2) == 3
    view(-35, 25); grid on; zlabel('w [m]');
else
    axis equal; set(gca, 'YDir', 'reverse');
end
xlabel('span x [m]'); ylabel('chord y [m]'); title(ttl);
axis tight;
end
