function d = dof_index(nodes, comps)
%DOF_INDEX  Global DOF numbers.  comps: 1=u 2=v 3=w 4=psi_x 5=psi_y.
%   d = dof_index(nodes, comps) returns a column vector for all
%   combinations of the given nodes and components.
nodes = nodes(:); comps = comps(:)';
d = 5*(nodes - 1) + comps;
d = d(:);
end
