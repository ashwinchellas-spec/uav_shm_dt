function out = navier_fsdt(sec, a, b, mode, varargin)
%NAVIER_FSDT  Closed-form Navier solution of FSDT for simply supported
%             cross-ply laminates (SS-1 conditions, Reddy 2004, Ch. 6-7).
%
%   out = navier_fsdt(sec, a, b, 'static', load, q0, nterms)
%       load = 'sin' (q = q0 sin(pi x/a) sin(pi y/b)) | 'uniform'
%       out.w(x,y), out.fields(x,y,z) handles; out.wc centre deflection
%   out = navier_fsdt(sec, a, b, 'modal', nmax)
%       out.omega (sorted, rad/s), out.mn (wave numbers of each)
%
%   Expansion (satisfies SS-1 exactly):
%     u = U cos(al x) sin(be y),  v = V sin(al x) cos(be y),
%     w = W sin(al x) sin(be y),  psi_x = X cos sin, psi_y = Y sin cos
%   with al = m pi / a, be = n pi / b.  Substitution into the five FSDT
%   equilibrium equations gives  ([S] - omega^2 [Mm]) {U V W X Y} = {0 0 q 0 0}
%     S11 = A11 al^2 + A66 be^2      S12 = (A12 + A66) al be
%     S14 = B11 al^2 + B66 be^2      S15 = (B12 + B66) al be
%     S22 = A66 al^2 + A22 be^2      S24 = S15
%     S25 = B66 al^2 + B22 be^2      S33 = A55 al^2 + A44 be^2
%     S34 = A55 al                   S35 = A44 be
%     S44 = D11 al^2 + D66 be^2 + A55
%     S45 = (D12 + D66) al be        S55 = D66 al^2 + D22 be^2 + A44
%   Valid only when A16 = A26 = B16 = B26 = D16 = D26 = A45 = 0 (cross-ply).
%   Note the shear indices: A55 <-> xz (our As(1,1)), A44 <-> yz (As(2,2)).

A = sec.A; B = sec.B; D = sec.D;
A55 = sec.As(1,1); A44 = sec.As(2,2);
Mm = [sec.I0 0 0 sec.I1 0; 0 sec.I0 0 0 sec.I1; 0 0 sec.I0 0 0;
      sec.I1 0 0 sec.I2 0; 0 sec.I1 0 0 sec.I2];
Smat = @(al, be) sym5(A, B, D, A55, A44, al, be);

switch lower(mode)
    case 'static'
        load = varargin{1}; q0 = varargin{2};
        if numel(varargin) > 2, nt = varargin{3}; else, nt = 99; end
        if strcmpi(load, 'sin'), ms = 1; else, ms = 1:2:nt; end
        C = zeros(5, numel(ms), numel(ms)); MN = zeros(2, numel(ms), numel(ms));
        for i = 1:numel(ms)
            for j = 1:numel(ms)
                m = ms(i); n = ms(j);
                if strcmpi(load, 'sin'), Q = q0; else, Q = 16*q0/(pi^2*m*n); end
                C(:,i,j) = Smat(m*pi/a, n*pi/b) \ [0;0;Q;0;0];
                MN(:,i,j) = [m; n];
            end
        end
        C = reshape(C, 5, []); MN = reshape(MN, 2, []);
        out.coef = C; out.mn = MN;
        out.w = @(x, y) navier_field(C, MN, a, b, x, y, 3);
        out.fields = @(x, y, z) navier_stress(sec, C, MN, a, b, x, y, z);
        out.wc = out.w(a/2, b/2);
    case 'modal'
        nmax = varargin{1};
        om = []; mn = [];
        for m = 1:nmax
            for n = 1:nmax
                S = Smat(m*pi/a, n*pi/b);
                lam = sort(real(eig(S, Mm)));
                om = [om; sqrt(lam(1))]; mn = [mn; m n];
            end
        end
        [out.omega, i] = sort(om); out.mn = mn(i,:);
    otherwise
        error('navier_fsdt: unknown mode');
end
end

function S = sym5(A, B, D, A55, A44, al, be)
S = zeros(5);
S(1,1) = A(1,1)*al^2 + A(3,3)*be^2;
S(1,2) = (A(1,2) + A(3,3))*al*be;
S(1,4) = B(1,1)*al^2 + B(3,3)*be^2;
S(1,5) = (B(1,2) + B(3,3))*al*be;
S(2,2) = A(3,3)*al^2 + A(2,2)*be^2;
S(2,4) = S(1,5);
S(2,5) = B(3,3)*al^2 + B(2,2)*be^2;
S(3,3) = A55*al^2 + A44*be^2;
S(3,4) = A55*al;
S(3,5) = A44*be;
S(4,4) = D(1,1)*al^2 + D(3,3)*be^2 + A55;
S(4,5) = (D(1,2) + D(3,3))*al*be;
S(5,5) = D(3,3)*al^2 + D(2,2)*be^2 + A44;
S = S + triu(S,1)';
end

function w = navier_field(C, MN, a, b, x, y, comp)
al = MN(1,:)*pi/a; be = MN(2,:)*pi/b;
w = zeros(size(x));
for k = 1:size(C,2)
    w = w + C(comp,k) * sin(al(k)*x) .* sin(be(k)*y);
end
end

function f = navier_stress(sec, C, MN, a, b, x, y, z)
% in-plane global stresses at (x,y,z) and shear resultants at (x,y)
al = MN(1,:)*pi/a; be = MN(2,:)*pi/b;
e0 = zeros(3,1); kap = zeros(3,1); Q = zeros(2,1);
for k = 1:size(C,2)
    U = C(1,k); V = C(2,k); W = C(3,k); X = C(4,k); Y = C(5,k);
    ss = sin(al(k)*x)*sin(be(k)*y); cc = cos(al(k)*x)*cos(be(k)*y);
    cs = cos(al(k)*x)*sin(be(k)*y); sc = sin(al(k)*x)*cos(be(k)*y);
    e0  = e0  + [-al(k)*U*ss; -be(k)*V*ss; (be(k)*U + al(k)*V)*cc];
    kap = kap + [-al(k)*X*ss; -be(k)*Y*ss; (be(k)*X + al(k)*Y)*cc];
    Q   = Q   + [sec.As(1,1)*(al(k)*W + X)*cs; sec.As(2,2)*(be(k)*W + Y)*sc];
end
kply = find(z >= sec.z(1:end-1) - 1e-15 & z <= sec.z(2:end) + 1e-15, 1, 'last');
f.sig = sec.Qb{kply} * (e0 + z*kap);
f.e0 = e0; f.kap = kap; f.Q = Q;
end
