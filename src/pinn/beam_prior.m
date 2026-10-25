function pr = beam_prior(cfg, n_ult, x_off, nr)
%BEAM_PRIOR  Reduced-order physics model of the wing (the "physics" half of
%   the hybrid physics-AI surrogate).
%
%   pr = beam_prior(cfg, n_ult, x_off, nr)
%   cfg   : wing configuration (sample_to_case output)
%   n_ult : ultimate load factor;  x_off : stress read-out station [m]
%   nr    : number of Ritz functions (default 6)
%
%   The plate wing is condensed to a tapered cantilever beam on an elastic
%   root using Classical Lamination Theory for the section:
%     bending stiffness  EI(x) = c(x) / d11,   d = D^-1 (laminate compliance:
%                        free anticlastic curvature, as in a wide beam)
%     shear stiffness    kGA(x) = c(x) As11  (energy-equivalent k, FSDT)
%     torsion            GJ(x) = 4 D66 c(x)   (Kirchhoff plate in torsion)
%     mass / length      m(x) = I0 c(x),   polar inertia  m c^2/12
%     root springs       K_theta = k* D11 / L * c_root  (as in the FE model)
%   1. f_1B : Rayleigh-Ritz, w = sum a_j (x/L)^j, j = 1..nr (j = 1 is the
%      rigid rotation allowed by the root spring), corrected for transverse
%      shear with Dunkerley's rule 1/f^2 = 1/f_b^2 + 1/f_s^2,
%      f_s = (1/4L) sqrt(kGA/m)  (first shear-beam mode).
%   2. f_1T : uniform-torsion estimate (1/4L) sqrt(GJ/(m c^2/12)) with the
%      mean chord, root torsion spring in series (Dunkerley).
%   3. w_tip: static Ritz solution under the net 1 g load of the FE model
%      (Schrenk lift W/2 minus wing weight) + shear deflection int V/kGA.
%   4. RF_TW: Tsai-Wu strength ratio of the face plies at x_off under n_ult
%      times the 1 g bending moment (CLT ply stresses, first-ply failure).
%   pr.y = [f1B f1T wtip RF];  pr.feat: dimensionless CLT groups (below).
%   It ignores sweep (bending-torsion coupling), chordwise bending,
%   2-D load diffusion near the root and the true pressure distribution:
%   that missing physics is what the neural network has to learn.
if nargin < 4, nr = 6; end
g = cfg.geom; L = g.span; c0 = g.c_root; lam = g.taper;
sec = laminate_section(wing_layup(cfg), cfg.lam.kshear);
d = inv([sec.A sec.B; sec.B sec.D]);
d11 = d(4,4);
xg = linspace(0, L, 401)';
c = c0*(1 - (1 - lam)*xg/L);
EI = c/d11; kGA = c*sec.As(1,1); m = sec.I0*c;
Kth = cfg.fe.kroot*sec.D(1,1)/L*c0;
% Ritz basis (x/L)^j and derivatives
j = 1:nr; s = xg/L;
Ph = s.^j; Ph2 = (j.*(j-1)).*s.^max(j-2,0)/L^2; Ph2(:, j < 2) = 0;
dPh0 = (j == 1)/L;
wq = trapz_w(xg);
K = Ph2'*(wq.*EI.*Ph2) + Kth*(dPh0'*dPh0);
M = Ph'*(wq.*m.*Ph);
lb = sort(real(eig((K + K')/2, (M + M')/2)));
fb = sqrt(lb(1))/(2*pi);
mb = sum(wq.*m)/L;
fs = 1/(4*L)*sqrt(mean(kGA)/mb);
pr.f1B = 1/sqrt(1/fb^2 + 1/fs^2);
% torsion
cm = mean(c); GJ = 4*sec.D(3,3)*cm; Ip = sec.I0*cm^3/12;
ft = 1/(4*L)*sqrt(GJ/Ip);
Kt = cfg.fe.kroot*sec.D(2,2)/L*c0;                    % root torsion spring (as FE)
ftr = 1/(2*pi)*sqrt(Kt/(Ip*L/3));                     % rigid rotation on the spring
pr.f1T = 1/sqrt(1/ft^2 + 1/ftr^2);
% static 1 g net load
W = cfg.ac.mass*cfg.g;
S = trapz(xg, c);
p = W/2*0.5*(c/S + 4/(pi*L)*sqrt(max(1 - s.^2, 0))) - cfg.g*m;
a = K \ (Ph'*(wq.*p));
V = flipud(cumtrapz(flipud(-xg), flipud(p)));      % V(x) = int_x^L p
pr.wtip = Ph(end,:)*a + trapz(xg, V./kGA);
% Tsai-Wu at x_off under the ultimate load
Mx = trapz(xg, p.*max(xg - x_off, 0));
cx = c0*(1 - (1 - lam)*x_off/L);
mx = n_ult*Mx/cx;                                     % moment per unit width
ek = d*[0;0;0; mx;0;0];
R = inf;
for k = 1:sec.nply
    if strcmp(sec.layup.mat(k).type, 'core'), continue; end
    for zz = sec.z(k:k+1)
        sm = sec.Qm{k}*(sec.Te{k}*(ek(1:3) + zz*ek(4:6)));
        fi = failure_indices(sm, sec.layup.mat(k));
        R = min(R, 1/fi.tsaiwu);
    end
end
pr.RF = R;
pr.y = [pr.f1B, pr.f1T, pr.wtip, pr.RF];
% dimensionless CLT groups the beam model ignores (inputs for the network):
% orthotropy, Poisson coupling, bending-twist coupling (D16, D26), shear
% flexibility of the core relative to bending (root fixity, taper and sweep
% are design inputs already)
Dm = sec.D;
pr.feat = [log(Dm(2,2)/Dm(1,1)), log(Dm(3,3)/Dm(1,1)), Dm(1,2)/Dm(1,1), ...
           Dm(1,3)/sqrt(Dm(1,1)*Dm(3,3)), Dm(2,3)/sqrt(Dm(2,2)*Dm(3,3)), ...
           log(Dm(1,1)/(sec.As(1,1)*L^2))];
pr.feat_names = {'ln D22/D11', 'ln D66/D11', 'D12/D11', 'D16/sqrt(D11 D66)', 'D26/sqrt(D22 D66)', 'ln D11/(kGA L^2)'};
end

function w = trapz_w(x)
% trapezoidal quadrature weights
dx = diff(x(:)); w = zeros(numel(x),1);
w(1:end-1) = w(1:end-1) + dx/2; w(2:end) = w(2:end) + dx/2;
end
