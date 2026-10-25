function out = modal_response(omega, zeta, Phi_s, gam, g, dt, resid_s)
%MODAL_RESPONSE  Exact-discretised modal superposition (mode-acceleration form).
%
%   out = modal_response(omega, zeta, Phi_s, gam, g, dt, resid_s)
%   Load F(t) = F0 g(t). Each mode obeys
%       q'' + 2 zeta w q' + w^2 q = gam_i g(t),   gam = Phi' F0
%   and is integrated EXACTLY for a piecewise-linear g (first-order hold):
%       x = [q; q'],  x_{k+1} = E x_k + G1 g_k + G2 (g_{k+1} - g_k)
%   with  expm([A B 0; 0 0 1/dt; 0 0 0] dt) = [E G1 G2; 0 1 1; 0 0 1]
%   (Franklin & Powell).  No period error, unconditionally stable.
%
%   Mode-acceleration correction: the static contribution of all modes
%   that were NOT kept is added exactly,
%       d(t) = Phi_s q(t) + resid_s g(t),
%       resid_s = (K^-1 F0 - sum_i phi_i gam_i / w_i^2) at the output DOFs.
%
%   omega, zeta : nm x 1   (zeta scalar allowed)
%   Phi_s       : ns x nm  mode shapes at the output DOFs (or any linear
%                 output, e.g. modal strains at sensor locations)
%   gam         : nm x 1   modal forces
%   g           : 1 x nt   load history samples, uniform step dt
%   resid_s     : ns x 1   residual static flexibility term (optional)
%   out.d (ns x nt) displacement-type output, out.a (ns x nt) acceleration
%   (the static correction carries no acceleration), out.q modal coords.
nm = numel(omega); nt = numel(g);
if isscalar(zeta), zeta = zeta*ones(nm,1); end
q = zeros(nm, nt); qd = q; qdd = q;
for i = 1:nm
    w = omega(i); z = zeta(i);
    A = [0 1; -w^2 -2*z*w]; B = [0; gam(i)];
    Mx = zeros(4); Mx(1:2,1:2) = A; Mx(1:2,3) = B; Mx(3,4) = 1/dt;
    Ex = expm(Mx*dt);
    E = Ex(1:2,1:2); G1 = Ex(1:2,3); G2 = Ex(1:2,4);
    x = [0; 0];
    for k = 1:nt-1
        x = E*x + G1*g(k) + G2*(g(k+1) - g(k));
        q(i,k+1) = x(1); qd(i,k+1) = x(2);
    end
    qdd(i,:) = gam(i)*g - 2*z*w*qd(i,:) - w^2*q(i,:);
end
out.d = Phi_s * q;
if nargin >= 7 && ~isempty(resid_s), out.d = out.d + resid_s(:) * g; end
out.a = Phi_s * qdd;
out.q = q;
end
