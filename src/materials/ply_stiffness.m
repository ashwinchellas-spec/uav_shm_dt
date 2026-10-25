function [Q, Qs] = ply_stiffness(mat)
%PLY_STIFFNESS  Plane-stress reduced stiffness of a lamina (material axes).
%
%   [Q, Qs] = ply_stiffness(mat)
%   Q  (3x3) relates [s11 s22 t12]' = Q [e11 e22 g12]'
%   Qs (2x2) relates [t13 t23]'     = Qs [g13 g23]'
%
%   Q11 = E1/(1-nu12 nu21), Q22 = E2/(1-nu12 nu21),
%   Q12 = nu12 E2/(1-nu12 nu21), Q66 = G12, nu21 = nu12 E2/E1.

nu21  = mat.nu12 * mat.E2 / mat.E1;
den   = 1 - mat.nu12 * nu21;
Q11   = mat.E1 / den;
Q22   = mat.E2 / den;
Q12   = mat.nu12 * mat.E2 / den;
Q66   = mat.G12;
Q  = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
Qs = [mat.G13 0; 0 mat.G23];
end
