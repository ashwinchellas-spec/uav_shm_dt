function [pfun, Lp] = lift_distribution(geom, Ltot, spanwise)
%LIFT_DISTRIBUTION  Aerodynamic pressure on a semi-wing (physics-based).
%
%   [pfun, Lp] = lift_distribution(geom, Ltot, spanwise)
%   geom     : wing geometry (span, c_root, taper, sweep_le)
%   Ltot     : lift carried by the semi-wing [N]
%   spanwise : 'elliptic' | 'planform' | 'schrenk' (default)
%
%   Spanwise lift per unit span L'(x), x in [0, b]:
%     elliptic : L'_e = (4 Ltot / (pi b)) sqrt(1 - (x/b)^2)
%     planform : L'_p = Ltot c(x) / S_semi
%     Schrenk  : L'   = (L'_e + L'_p) / 2   (Schrenk, NACA TM-948, 1940)
%
%   Chordwise loading (thin-airfoil theory, flat plate at incidence):
%     Delta_p(s) ~ sqrt((1-s)/s),  s = (y - y_le)/c,
%   normalised  f(s) = (2/pi) sqrt((1-s)/s),  int_0^1 f ds = 1,
%   whose centroid is exactly the quarter chord (aerodynamic centre), so
%   the load also twists the wing about its elastic axis.
%
%   p(x,y) = L'(x) f(s) / c(x)          [Pa]
%   Lp     : handle for L'(x)

if nargin < 3, spanwise = 'schrenk'; end
if ~isfield(geom, 'taper'),    geom.taper = 1;    end
if ~isfield(geom, 'sweep_le'), geom.sweep_le = 0; end
b  = geom.span; cr = geom.c_root; lam = geom.taper;
S  = b * cr * (1 + lam) / 2;
cf = @(x) cr * (1 - (1 - lam) * x / b);
Le = @(x) 4*Ltot/(pi*b) * sqrt(max(0, 1 - (x/b).^2));
Lpl = @(x) Ltot * cf(x) / S;
switch lower(spanwise)
    case 'elliptic', Lp = Le;
    case 'planform', Lp = Lpl;
    case 'schrenk',  Lp = @(x) 0.5*(Le(x) + Lpl(x));
    otherwise, error('lift_distribution: unknown option');
end
fch  = @(s) (2/pi) * sqrt(max(0, 1 - s) ./ max(s, 1e-12));
sfun = @(x, y) (y - x*tand(geom.sweep_le)) ./ cf(x);
pfun = @(x, y) Lp(x) .* fch(sfun(x, y)) ./ cf(x);
end
