function yt = naca4_thickness(s, tc)
%NACA4_THICKNESS  Half-thickness of a symmetric NACA 4-digit airfoil (per unit chord).
%   yt = 5 tc (0.2969 sqrt(s) - 0.1260 s - 0.3516 s^2 + 0.2843 s^3 - 0.1015 s^4)
yt = 5*tc*(0.2969*sqrt(s) - 0.1260*s - 0.3516*s.^2 + 0.2843*s.^3 - 0.1015*s.^4);
end
