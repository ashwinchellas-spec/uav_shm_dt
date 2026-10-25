function [label, frac] = classify_mode(model, phi)
%CLASSIFY_MODE  Label a mode shape by its dominant deformation.
%   [label, frac] = classify_mode(model, phi)
%   Uses the w-field along the span and chord:
%     bending  - w nearly constant along every chord line
%     torsion  - w varies linearly along the chord (rotation about x)
%     chordwise/plate - higher chord-wise curvature
%     in-plane - in-plane (u, v) kinetic energy dominates
%   The number of nodal lines along the span gives the order (1B, 2B, 1T ...).
%   frac = [bending torsion chordwise in-plane] fractions of the w-energy.
mesh = model.mesh; nx = mesh.nx; ny = mesh.ny;
w = reshape(phi(3:5:end), nx+1, ny+1);           % (span, chord)
uv = phi([1:5:end, 2:5:end]);
eta = linspace(-1, 1, ny+1);
P = [ones(ny+1,1) eta(:) (3*eta(:).^2 - 1)/2];   % Legendre basis along chord
c = w / P';                                       % least squares per span station
res = w - c * P';
e = [sum(c(:,1).^2)*2, sum(c(:,2).^2)*2/3, sum(c(:,3).^2)*2/5 + sum(res(:).^2)*2/(ny+1)];
wtot = sum(w(:).^2) + eps;
inpl = sum(uv.^2) / (sum(uv.^2) + wtot);
frac = [e / sum(e) * (1 - inpl), inpl];
[~, k] = max(frac);
switch k
    case 1, s = c(:,1); base = 'B';
    case 2, s = c(:,2); base = 'T';
    case 3, s = c(:,3); base = 'C';
    case 4
        label = 'in-plane'; return;
end
s = s(2:end);                                     % skip clamped root
nodal = sum(abs(diff(sign(s(abs(s) > 1e-3*max(abs(s)))))) > 0);
label = sprintf('%d%s', nodal + 1, base);
end
