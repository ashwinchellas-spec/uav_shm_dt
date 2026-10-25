function sc = shear_correction(sec)
%SHEAR_CORRECTION  Energy-equivalent FSDT shear correction factors.
%
%   sc = shear_correction(sec)  (sec from laminate_section)
%
%   For cylindrical bending in direction x the equilibrium equation
%   d(tau_xz)/dz = -d(sigma_x)/dx gives the exact layer-wise shear stress
%
%       tau_xz(z) = -(Q_x / D_x) g(z),   g(z) = int_{-h/2}^{z} E(s)(s - z_n) ds
%
%   with E = Qb11(z), z_n the neutral axis and D_x = int E (z - z_n)^2 dz.
%   Equating the shear strain energy of this distribution with the FSDT
%   (constant shear strain) energy  Q_x^2 / (2 k int G dz)  gives
%
%       k_x = D_x^2 / ( int G dz  *  int g(z)^2 / G(z) dz ),   G = Qsb55.
%
%   For a homogeneous plate this reproduces k = 5/6 exactly.  For sandwich
%   panels it is far below 5/6 and k = 5/6 grossly over-predicts the
%   transverse shear stiffness. (Whitney 1973; Vlachoutsis 1992.)
%
%   sc.kx, sc.ky : factors
%   sc.tprof.z   : sample z-points, 5 per ply (interfaces included)
%   sc.tprof.fx  : tau_xz(z)/Q_x   (integrates to 1 over the thickness)
%   sc.tprof.fy  : tau_yz(z)/Q_y

N = sec.nply; z = sec.z;
Ex = zeros(1,N); Ey = Ex; Gx = Ex; Gy = Ex;
for k = 1:N
    Ex(k) = sec.Qb{k}(1,1);  Ey(k) = sec.Qb{k}(2,2);
    Gx(k) = sec.Qsb{k}(1,1); Gy(k) = sec.Qsb{k}(2,2);
end
[sc.kx, zx, fx] = one_dir(Ex, Gx, z);
[sc.ky, ~,  fy] = one_dir(Ey, Gy, z);
sc.tprof.z = zx; sc.tprof.fx = fx; sc.tprof.fy = fy;
end

function [k, zs, f] = one_dir(E, G, z)
N  = numel(E);
t  = diff(z);
zn = sum(E .* (z(2:end).^2 - z(1:end-1).^2) / 2) / sum(E .* t);
Dn = sum(E .* ((z(2:end)-zn).^3 - (z(1:end-1)-zn).^3) / 3);
% 3-point Gauss is exact for the quartic g^2 within each ply
gp = [-sqrt(3/5) 0 sqrt(3/5)]; gw = [5 8 5]/9;
gfun = @(zz, gb, Ek, zb) gb + Ek * ((zz.^2 - zb^2)/2 - zn*(zz - zb));
gb = 0; Ig2 = 0;
zs = zeros(1, 5*N); f = zs;
for kk = 1:N
    zb = z(kk); zt = z(kk+1);
    zq = (zb+zt)/2 + gp*(zt-zb)/2;
    Ig2 = Ig2 + sum(gw .* gfun(zq, gb, E(kk), zb).^2) * (zt-zb)/2 / G(kk);
    zl = linspace(zb, zt, 5);
    zs(5*kk-4:5*kk) = zl;
    f(5*kk-4:5*kk)  = -gfun(zl, gb, E(kk), zb) / Dn;
    gb = gfun(zt, gb, E(kk), zb);
end
k = Dn^2 / (sum(G .* t) * Ig2);
end
