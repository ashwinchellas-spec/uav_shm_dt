function fi = failure_indices(s, mat)
%FAILURE_INDICES  Lamina failure indices from material-axis stresses.
%
%   fi = failure_indices(s, mat)
%   s   : 3 x n stresses [s11; s22; t12] in material axes [Pa]
%   mat : material struct with strengths Xt, Xc, Yt, Yc, S12, S23
%
%   fi.tsaiwu   : Tsai-Wu index expressed as 1/R, where R is the strength
%                 ratio (load multiplier to failure).  >= 1 means failure.
%                 F12* = -1/2 sqrt(F11 F22)  (Tsai & Hahn)
%   fi.maxstress: maximum-stress index
%   fi.hashin   : 4 x n  [fibre tension; fibre compression;
%                         matrix tension; matrix compression]  (Hashin 1980)
%   fi.hashin_max, fi.mode (1..4) : governing Hashin mode

s1 = s(1,:); s2 = s(2,:); t12 = s(3,:);
Xt = mat.Xt; Xc = mat.Xc; Yt = mat.Yt; Yc = mat.Yc; S = mat.S12; ST = mat.S23;

F1 = 1/Xt - 1/Xc; F2 = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc); F22 = 1/(Yt*Yc); F66 = 1/S^2; F12 = -0.5*sqrt(F11*F22);
a = F11*s1.^2 + F22*s2.^2 + F66*t12.^2 + 2*F12*s1.*s2;   % quadratic part
b = F1*s1 + F2*s2;                                         % linear part
% a R^2 + b R - 1 = 0  ->  positive root R; index = 1/R
R = zeros(size(a));
q = a > 0;
R(q)  = (-b(q) + sqrt(b(q).^2 + 4*a(q))) ./ (2*a(q));
R(~q & b > 0) = 1 ./ b(~q & b > 0);
R(~q & b <= 0) = inf;
fi.tsaiwu = 1 ./ R;

sx = s1; sx(s1 >= 0) = s1(s1 >= 0)/Xt; sx(s1 < 0) = -s1(s1 < 0)/Xc;
sy = s2; sy(s2 >= 0) = s2(s2 >= 0)/Yt; sy(s2 < 0) = -s2(s2 < 0)/Yc;
fi.maxstress = max([sx; sy; abs(t12)/S], [], 1);

n = numel(s1); H = zeros(4, n);
ft = s1 >= 0; mt = s2 >= 0;
H(1, ft)  = (s1(ft)/Xt).^2 + (t12(ft)/S).^2;
H(2, ~ft) = (s1(~ft)/Xc).^2;
H(3, mt)  = (s2(mt)/Yt).^2 + (t12(mt)/S).^2;
H(4, ~mt) = (s2(~mt)/(2*ST)).^2 + ((Yc/(2*ST))^2 - 1)*s2(~mt)/Yc + (t12(~mt)/S).^2;
H = sqrt(max(H, 0));        % stress-like (linear) scaling, 1 = failure
fi.hashin = H;
[fi.hashin_max, fi.mode] = max(H, [], 1);
end
