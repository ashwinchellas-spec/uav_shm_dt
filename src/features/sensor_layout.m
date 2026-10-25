function sens = sensor_layout(ref, n_acc, fbg_lines, n_grat)
%SENSOR_LAYOUT  Realistic SHM sensor network of the digital twin.
%
%   sens = sensor_layout(ref, n_acc, fbg_lines, n_grat)
%   1. FBG strain sensing (bottom skin, span-wise strain eps_x):
%      optical fibres along the span at chord fractions fbg_lines (e.g.
%      [0.2 0.5 0.8]), each with n_grat equally spaced Bragg gratings from
%      x = 0.05 L to 0.95 L.  Static (tip-test) and dynamic (modal strain)
%      readings.
%   2. Accelerometers (z): n_acc locations chosen by Effective Independence
%      on the 5 flexural reference modes.  Candidates form a practical
%      mounting grid - span stations every 0.1 L on the leading edge,
%      mid-chord and trailing edge (30 points).  Unconstrained EfI on all
%      nodes clusters sensors on ADJACENT nodes (it assumes independent
%      noise, so near-duplicate information looks attractive); the grid
%      enforces >= 100 mm spacing, as in practice.
%   3. One temperature sensor.
mesh = ref.mesh; L = mesh.geom.span; c = mesh.geom.c_root;
xc = elem_centroids(mesh);
xs = linspace(0.05, 0.95, n_grat)*L;
fe = zeros(numel(fbg_lines)*n_grat, 1); pos = zeros(numel(fe), 2); k = 0;
for a = 1:numel(fbg_lines)
    for b = 1:n_grat
        k = k + 1;
        [~, fe(k)] = min((xc(:,1) - xs(b)).^2 + (xc(:,2) - fbg_lines(a)*c).^2);
        pos(k,:) = xc(fe(k),:);
    end
end
sens.fbg_elem = fe; sens.fbg_pos = pos; sens.fbg_lines = fbg_lines; sens.n_grat = n_grat;
cand = [];
for xs_ = (0.1:0.1:1.0)*L
    for yc_ = [0 0.5 1]*c
        [~, nd] = min((mesh.nodes(:,1) - xs_).^2 + (mesh.nodes(:,2) - yc_).^2);
        cand(end+1,1) = nd; %#ok<AGROW>
    end
end
sens.acc_cand = cand;
Phi = double(ref.o.Wmodes(cand, ref.flex));
[sens.acc_node, sens.efi_hist] = efi_placement(Phi, cand, n_acc);
sens.acc_pos = mesh.nodes(sens.acc_node, :);
P = double(ref.o.Wmodes(sens.acc_node, ref.flex));
sens.fisher_cond = cond(P'*P);
% MAC between the reference modes seen by the accelerometers: should be ~0
% off-diagonal (modes distinguishable with this sensor set)
M = (P'*P).^2 ./ (diag(P'*P)*diag(P'*P)');
sens.automac_offdiag = max(M(~eye(size(M))));
end
