function X = design_inputs(P, spec)
%DESIGN_INPUTS  Encode design-tier variables as network inputs.
%   Ply angles are 180 deg periodic -> [cos 2theta, sin 2theta] (a 0 and a
%   180 deg ply are the same ply); face material -> one-hot; log-scaled
%   variables (core grade, root fixity) -> log10; all others linear.
names = {spec.vars.name};
v = @(nm) P(:, strcmp(names, nm));
X = [];
for k = 1:4
    th = v(sprintf('th%d', k));
    X = [X, cosd(2*th), sind(2*th)]; %#ok<AGROW>
end
im = min(numel(spec.face_mats), floor(v('face_mat')*numel(spec.face_mats)) + 1);
X = [X, v('tply_mm'), v('tcore_mm'), log10(v('core_grade')), double(im == 1:numel(spec.face_mats)), ...
     v('taper'), v('sweep'), log10(v('kroot'))];
end
