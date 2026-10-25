function xc = elem_centroids(mesh)
%ELEM_CENTROIDS  ne x 2 element centroid coordinates (average of the 4 nodes).
xc = zeros(mesh.ne, 2);
for a = 1:4, xc = xc + mesh.nodes(mesh.elems(:,a), :) / 4; end
end
