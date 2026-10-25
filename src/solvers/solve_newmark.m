function out = solve_newmark(model, F0, gfun, t, opts)
%SOLVE_NEWMARK  Implicit Newmark-beta transient analysis.
%
%   out = solve_newmark(model, F0, gfun, t, opts)
%   Solves  M a + C v + K d = F0 * g(t)   on the free DOFs.
%   F0    : spatial load vector (ndof x 1)
%   gfun  : handle g(t) (load time history), vectorised
%   t     : uniform time vector
%   opts.alpha, opts.beta_r : Rayleigh damping  C = alpha M + beta_r K
%   opts.gamma, opts.beta   : Newmark parameters (default 1/2, 1/4:
%                             average acceleration, unconditionally stable,
%                             no numerical damping)
%   opts.out_dofs           : DOFs to store (default all)
%   opts.d0, opts.v0        : initial conditions (default zero)
%
%   out.t, out.d, out.v, out.a (stored DOFs x nt), out.dofs, out.energy
%   (kinetic + strain energy, for verification).

if nargin < 5, opts = struct(); end
def = struct('alpha',0,'beta_r',0,'gamma',0.5,'beta',0.25, ...
             'out_dofs',(1:model.ndof)','d0',[],'v0',[]);
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end, end

f = model.free;
K = model.K(f,f); M = model.M(f,f);
C = opts.alpha * M + opts.beta_r * K;
nf = numel(f); nt = numel(t); dt = t(2) - t(1);
gm = opts.gamma; bt = opts.beta;

d = zeros(nf,1); v = zeros(nf,1);
if ~isempty(opts.d0), d = opts.d0(f); end
if ~isempty(opts.v0), v = opts.v0(f); end
Ff = F0(f); g = gfun(t);
a = M \ (Ff*g(1) - C*v - K*d);

a0 = 1/(bt*dt^2); a1 = gm/(bt*dt); a2 = 1/(bt*dt); a3 = 1/(2*bt) - 1;
a4 = gm/bt - 1;   a5 = dt/2*(gm/bt - 2);
Keff = K + a0*M + a1*C;
[R, p, P] = chol(Keff);           % R'R = P' Keff P  (sparse, reordered)
if p ~= 0, error('solve_newmark: effective stiffness not SPD'); end

[~, loc] = ismember(opts.out_dofs, f);
keep = loc > 0;
out.dofs = opts.out_dofs(:);
out.d = zeros(numel(out.dofs), nt); out.v = out.d; out.a = out.d;
out.energy = zeros(1, nt);
out.d(keep,1) = d(loc(keep));
out.v(keep,1) = v(loc(keep));
out.a(keep,1) = a(loc(keep));
out.energy(1) = 0.5*(v'*M*v) + 0.5*(d'*K*d);
for n = 1:nt-1
    rhs = Ff*g(n+1) + M*(a0*d + a2*v + a3*a) + C*(a1*d + a4*v + a5*a);
    dn  = P * (R \ (R' \ (P' * rhs)));
    an  = a0*(dn - d) - a2*v - a3*a;
    vn  = v + dt*((1-gm)*a + gm*an);
    d = dn; v = vn; a = an;
    out.d(keep,n+1) = d(loc(keep));
    out.v(keep,n+1) = v(loc(keep));
    out.a(keep,n+1) = a(loc(keep));
    out.energy(n+1) = 0.5*(v'*M*v) + 0.5*(d'*K*d);
end
out.t = t;
end
