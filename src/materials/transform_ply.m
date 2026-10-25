function [Qb, Qsb, Te] = transform_ply(Q, Qs, theta_deg)
%TRANSFORM_PLY  Rotate lamina stiffness from material (1,2) to global (x,y).
%
%   [Qb, Qsb, Te] = transform_ply(Q, Qs, theta_deg)
%   theta is the angle from the global x-axis (span) to the fibre axis 1,
%   positive counter-clockwise about z.
%
%   Qb  = T^-1 Q T^-T   (with engineering shear strains)  -> Q-bar
%   Qsb : transverse shear [t_xz t_yz]' = Qsb [g_xz g_yz]'
%   Te  : strain transformation, e_material = Te * e_global
%         (e = [ex ey gxy]'), used in stress recovery.

c = cosd(theta_deg); s = sind(theta_deg);

% Strain transformation (engineering shear), global -> material
Te = [ c^2   s^2    c*s;
       s^2   c^2   -c*s;
      -2*c*s 2*c*s  c^2 - s^2];
Qb = Te' * Q * Te;

% Transverse shear: [g13 g23]' = R [g_xz g_yz]',  R = [c s; -s c]
R   = [c s; -s c];
Qsb = R' * Qs * R;
end
