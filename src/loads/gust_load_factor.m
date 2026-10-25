function [dn, info] = gust_load_factor(ac, Ude)
%GUST_LOAD_FACTOR  Pratt-Walker gust load factor increment (CS-23 / FAR-23 style).
%
%   [dn, info] = gust_load_factor(ac, Ude)
%   ac.mass [kg], ac.S [m^2] (full wing area), ac.b_full [m], ac.cbar [m],
%   ac.V [m/s], ac.rho [kg/m^3]; Ude derived gust velocity [m/s].
%
%     a   = 2 pi AR / (2 + sqrt(AR^2 + 4))       (Helmbold lift slope)
%     mu  = 2 (W/S) / (rho cbar a g)            (mass ratio)
%     Kg  = 0.88 mu / (5.3 + mu)                (gust alleviation)
%     dn  = rho Ude V a Kg / (2 W/S)
g = 9.81;
W  = ac.mass * g;
AR = ac.b_full^2 / ac.S;
a  = 2*pi*AR / (2 + sqrt(AR^2 + 4));
WS = W / ac.S;
mu = 2*WS / (ac.rho * ac.cbar * a * g);
Kg = 0.88*mu / (5.3 + mu);
dn = ac.rho * Ude * ac.V * a * Kg / (2*WS);
info = struct('a', a, 'mu', mu, 'Kg', Kg, 'AR', AR, 'WS', WS);
end
