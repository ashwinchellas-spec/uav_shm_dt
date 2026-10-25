function fac = root_loosening_factor(mesh, dm)
%ROOT_LOOSENING_FACTOR  Root-spring multipliers for a loosened root attachment.
%   fac = root_loosening_factor(mesh, dm)   (one value per root node, sorted LE->TE)
%   The springs of the root nodes inside the chord band
%       |y - dm.center(2)| <= dm.radius(2)
%   are reduced to (1 - dm.severity) (bolt / adhesive loosening).
root = mesh.sets.root;
y = sort(mesh.nodes(root, 2));
fac = ones(numel(root), 1);
fac(abs(y - dm.center(2)) <= dm.radius(2) + 1e-12) = 1 - dm.severity;
end
