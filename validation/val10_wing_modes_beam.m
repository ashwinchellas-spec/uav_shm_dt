function R = val10_wing_modes_beam()
%VAL10_WING_MODES_BEAM  Independent check of the wing mode shapes & frequencies.
%
%   The 2-D FSDT plate modes of the baseline wing are compared with two
%   independent 1-D theories that share no code with the plate solver:
%
%   (a) BENDING (1B..5B): Timoshenko cantilever beam (own 2-node FE,
%       200 elements, shear + rotary inertia) with section properties
%       taken from the laminate:
%         EI   = c / d11,  d = D^-1  (free edges: anticlastic curvature allowed)
%         kGA  = c As11,   rhoA = I0 c,   rhoI = I2 c
%       Mode shapes compared with the Modal Assurance Criterion
%         MAC = (phi_p' phi_b)^2 / ((phi_p' phi_p)(phi_b' phi_b))
%       using the chord-averaged plate deflection along the span.
%
%   (b) TORSION (1T): Reissner (shear-deformable) plate-strip torsion
%         GJ = 4 D66 c [1 - tanh(lam c/2)/(lam c/2)],  lam = sqrt(As11/D66)
%         Ip = I0 c^3/12 + I2 c,   f_T1 = (1/4L) sqrt(GJ/Ip)
%       Kirchhoff plate theory (GJ = 4 D66 c) is printed for contrast.
%       Mode shape: twist theta(x) = (w_LE - w_TE)/c vs sin(pi x / 2L).
cfg = wing_baseline_config();
model = build_wing_model(cfg);
[fr, Phi] = solve_modal(model, 10);
labels = cell(10,1);
for k = 1:10, labels{k} = classify_mode(model, Phi(:,k)); end
sec = model.sections(1); c = cfg.geom.c_root; L = cfg.geom.span; mesh = model.mesh;
nx = mesh.nx; ny = mesh.ny;
xs = mesh.nodes(1:nx+1, 1);

% ---- (a) Timoshenko beam FE ------------------------------------------
dinv = inv(sec.D);
EI = c / dinv(1,1); kGA = c * sec.As(1,1); rA = sec.I0*c; rI = sec.I2*c;
ne = 200; le = L/ne; nd = 2*(ne+1);
Kb = zeros(nd); Mb = zeros(nd);
for e = 1:ne
    i = 2*e-1:2*e+2;                          % [w1 th1 w2 th2]
    kb = EI/le*[0 0 0 0; 0 1 0 -1; 0 0 0 0; 0 -1 0 1];
    Bs = [-1/le -0.5 1/le -0.5];              % gamma = w' - theta, 1-pt (no locking)
    ks = kGA*le*(Bs'*Bs);
    m1 = le/6*[2 1; 1 2];
    me = zeros(4); me([1 3],[1 3]) = rA*m1; me([2 4],[2 4]) = rI*m1;
    Kb(i,i) = Kb(i,i) + kb + ks; Mb(i,i) = Mb(i,i) + me;
end
fb = 3:nd;                                    % clamp w, theta at x = 0
[V, Lm] = eig(Kb(fb,fb), Mb(fb,fb));
[lam, o] = sort(diag(Lm)); V = V(:,o);
fbeam = sqrt(lam)/(2*pi);
xb = linspace(0, L, ne+1)';

R = [];
ib = find(~cellfun(@isempty, regexp(labels, 'B$')));
for j = 1:min(5, numel(ib))
    k = ib(j);
    wp = mean(reshape(Phi(3:5:end,k), nx+1, ny+1), 2);   % chord-averaged
    wb = V(1:2:end, j); wb = [0; wb];
    wbi = interp1(xb, wb, xs);
    mac = (wp'*wbi)^2 / ((wp'*wp)*(wbi'*wbi));
    R = [R, vres('V10 wing modes', sprintf('%s frequency, plate vs Timoshenko beam [Hz]', labels{k}), ...
         fr(k), fbeam(j), 3e-2), ...
         vres('V10 wing modes', sprintf('%s mode-shape MAC vs beam', labels{k}), mac, 1, 1e-2)]; %#ok<AGROW>
end

% ---- (b) Torsion -------------------------------------------------------
lamb = sqrt(sec.As(1,1)/sec.D(3,3));
phiR = 1 - tanh(lamb*c/2)/(lamb*c/2);
GJ = 4*sec.D(3,3)*c*phiR; GJk = 4*sec.D(3,3)*c;
Ip = sec.I0*c^3/12 + sec.I2*c;
fT1 = 1/(4*L)*sqrt(GJ/Ip); fT1k = 1/(4*L)*sqrt(GJk/Ip);
kt = find(strcmp(labels, '1T'), 1);
w2 = reshape(Phi(3:5:end,kt), nx+1, ny+1);
th = (w2(:,1) - w2(:,end))/c;
tb = sin(pi*xs/(2*L));
mac = (th'*tb)^2/((th'*th)*(tb'*tb));
R = [R, vres('V10 wing modes', '1T frequency, plate vs Reissner strip torsion [Hz]', fr(kt), fT1, 3e-2), ...
        vres('V10 wing modes', '1T mode-shape MAC vs sin(pi x/2L)', mac, 1, 1e-2), ...
        vres('V10 wing modes', sprintf('1T Kirchhoff torsion (no shear) [Hz] (Reissner factor %.3f)', phiR), ...
        fr(kt), fT1k, inf)];

% figure: plate vs beam mode shapes
fig = new_fig(1100, 420);
subplot(1,2,1); hold on; cl = lines(5);
for j = 1:min(4, numel(ib))
    k = ib(j);
    wp = mean(reshape(Phi(3:5:end,k), nx+1, ny+1), 2); wp = wp/max(abs(wp));
    wb = [0; V(1:2:end, j)]; wb = wb/max(abs(wb)) * sign(wb(end)) * sign(wp(end));
    plot(xs, wp, 'o', 'Color', cl(j,:), 'MarkerSize', 4);
    plot(xb, wb, '-', 'Color', cl(j,:), 'LineWidth', 1.3);
end
grid on; xlabel('span x [m]'); ylabel('normalised w');
title('Bending modes: plate FE (o) vs Timoshenko beam (-)');
subplot(1,2,2);
plot(xs, th/max(abs(th)), 'bo', xs, tb, 'b-', 'LineWidth', 1.3); grid on;
xlabel('span x [m]'); ylabel('normalised twist');
title(sprintf('1T: plate %.1f Hz, Reissner %.1f Hz, Kirchhoff %.1f Hz', fr(kt), fT1, fT1k));
save_fig(fig, 'V10_wing_modes_vs_beam');
end
