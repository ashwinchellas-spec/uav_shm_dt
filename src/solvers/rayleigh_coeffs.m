function [alpha, beta] = rayleigh_coeffs(w1, w2, z1, z2)
%RAYLEIGH_COEFFS  C = alpha M + beta K giving damping ratios z1, z2 at w1, w2 [rad/s].
%   zeta(w) = alpha/(2w) + beta w/2
if nargin < 4, z2 = z1; end
x = [1/(2*w1) w1/2; 1/(2*w2) w2/2] \ [z1; z2];
alpha = x(1); beta = x(2);
end
