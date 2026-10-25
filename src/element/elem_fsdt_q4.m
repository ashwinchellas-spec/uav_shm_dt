function [Ke, Me, rec] = elem_fsdt_q4(xy, sec, shear)
%ELEM_FSDT_Q4  4-node FSDT (Reissner-Mindlin) laminated plate element.
%
%   [Ke, Me, rec] = elem_fsdt_q4(xy, sec, shear)
%   xy    : 4x2 nodal coordinates (counter-clockwise)
%   sec   : section from laminate_section (A,B,D,As,I0,I1,I2)
%   shear : transverse-shear treatment
%     'hr'    (default) Hellinger-Reissner mixed formulation with
%             - an independently assumed covariant shear STRAIN field
%               gb_xi(eta)  = (1-eta)/2 g_xi(0,-1) + (1+eta)/2 g_xi(0,1)
%               gb_eta(xi)  = (1-xi)/2  g_eta(-1,0) + (1+xi)/2 g_eta(1,0)
%               (tied at the edge mid-points, where the displacement-based
%               strains are exact; MITC idea of Bathe & Dvorkin 1985)
%             - an independently assumed shear FORCE field, Piola mapped
%               Q = J' Qh / det J,   Qh = [b1 + b2 eta; b3 + b4 xi]
%             The HR functional  int( Q'gb - 1/2 Q' As^-1 Q ) dA  is then
%             stationary in beta for  beta = H^-1 G d  with
%               H = int P' As^-1 P dA,   G = int Ph' gb_cov dxi deta
%             (the Piola map makes G an exact parent-domain integral), so
%               Ks = G' H^-1 G.
%             Locking-free, passes the distorted patch test, exactly 6
%             rigid-body modes, and gives the shear forces directly.
%     'mitc4' assumed covariant shear strain only: Ks = int Bb_s' As Bb_s dA
%     'sri'   selective reduced integration (1-point shear; 2 spurious modes)
%     'full'  full 2x2 integration (locks for thin plates; reference only)
%
%   DOFs per node: [u v w psi_x psi_y]  ->  20 element DOFs.
%   Kinematics (FSDT):  u(x,y,z) = u0 + z psi_x,  v = v0 + z psi_y,  w = w0
%     membrane  e0    = [u,x ; v,y ; u,y + v,x]
%     bending   kappa = [psi_x,x ; psi_y,y ; psi_x,y + psi_y,x]
%     shear     gamma = [w,x + psi_x ; w,y + psi_y]
%   Constitutive: [N;M] = [A B; B D][e0; kappa],  Q = As gamma.
%
%   rec (element-centre recovery):  e0 = rec.Bm*de, kappa = rec.Bb*de,
%   Q = rec.Rq*de,  rec.area.

if nargin < 3 || isempty(shear), shear = 'hr'; end
shear = lower(shear);

C  = [sec.A sec.B; sec.B sec.D];
As = sec.As;
mI = [sec.I0 0 0 sec.I1 0;
      0 sec.I0 0 0 sec.I1;
      0 0 sec.I0 0 0;
      sec.I1 0 0 sec.I2 0;
      0 sec.I1 0 0 sec.I2];

% covariant shear-strain rows at the four tying points (A,C: g_xi; D,B: g_eta)
if any(strcmp(shear, {'hr', 'mitc4'}))
    [~, ~, ~, BsA, ~, ~, JA] = bmats(xy, 0, -1);
    [~, ~, ~, BsC, ~, ~, JC] = bmats(xy, 0,  1);
    [~, ~, ~, BsD, ~, ~, JD] = bmats(xy, -1, 0);
    [~, ~, ~, BsB, ~, ~, JB] = bmats(xy,  1, 0);
    gA = JA(1,:)*BsA; gC = JC(1,:)*BsC;     % g_xi  = [x,xi y,xi] . gamma
    gD = JD(2,:)*BsD; gB = JB(2,:)*BsB;     % g_eta = [x,eta y,eta] . gamma
end

[gp, gw] = gauss_1d(2);
Ke = zeros(20); Me = zeros(20); Ks = zeros(20);
H = zeros(4); G = zeros(4, 20);
area = 0;
for i = 1:2
    for j = 1:2
        xi = gp(i); eta = gp(j); w = gw(i)*gw(j);
        [N, Bm, Bb, Bs, Nm, dJ, J] = bmats(xy, xi, eta); %#ok<ASGLU>
        dA = dJ * w; area = area + dA;
        Bmb = [Bm; Bb];
        Ke = Ke + Bmb' * C * Bmb * dA;
        Me = Me + Nm' * mI * Nm * dA;
        switch shear
            case 'full'
                Ks = Ks + Bs' * As * Bs * dA;
            case {'hr', 'mitc4'}
                Bcov = [(1-eta)/2*gA + (1+eta)/2*gC;
                        (1-xi)/2*gD  + (1+xi)/2*gB];
                if strcmp(shear, 'mitc4')
                    Bsb = J \ Bcov;
                    Ks = Ks + Bsb' * As * Bsb * dA;
                else
                    Ph = [1 eta 0 0; 0 0 1 xi];
                    P  = J' * Ph / dJ;
                    H  = H + P' * (As \ P) * dA;
                    G  = G + Ph' * Bcov * w;
                end
        end
    end
end

[~, Bm0, Bb0, Bs0, ~, dJ0, J0] = bmats(xy, 0, 0);
switch shear
    case 'sri'
        Ks = Bs0' * As * Bs0 * (4 * dJ0);
        Rq = As * Bs0;
    case 'full'
        Rq = As * Bs0;
    case 'mitc4'
        Bcov0 = [(gA + gC)/2; (gD + gB)/2];
        Rq = As * (J0 \ Bcov0);
    case 'hr'
        Sb = H \ G;                              % beta = Sb * de
        Ks = G' * Sb;
        Rq = (J0' * [1 0 0 0; 0 0 1 0] / dJ0) * Sb;   % Q at the centre
    otherwise
        error('elem_fsdt_q4: unknown shear option "%s"', shear);
end
Ke = Ke + Ks;
Ke = (Ke + Ke') / 2;

if nargout > 2
    rec.Bm = Bm0; rec.Bb = Bb0; rec.Rq = Rq; rec.area = area; rec.Bs = Bs0;
end
end

% -------------------------------------------------------------------------
function [N, Bm, Bb, Bs, Nm, dJ, J] = bmats(xy, xi, eta)
[N, dNdxi] = shape_q4(xi, eta);
J  = dNdxi * xy;                     % [x,xi y,xi; x,eta y,eta]
dJ = det(J);
if dJ <= 0, error('elem_fsdt_q4: non-positive Jacobian (check node order)'); end
dN = J \ dNdxi;                      % [dN/dx; dN/dy]
Bm = zeros(3,20); Bb = zeros(3,20); Bs = zeros(2,20); Nm = zeros(5,20);
for a = 1:4
    c = 5*(a-1);
    Bm(1, c+1) = dN(1,a);
    Bm(2, c+2) = dN(2,a);
    Bm(3, c+1) = dN(2,a);  Bm(3, c+2) = dN(1,a);
    Bb(1, c+4) = dN(1,a);
    Bb(2, c+5) = dN(2,a);
    Bb(3, c+4) = dN(2,a);  Bb(3, c+5) = dN(1,a);
    Bs(1, c+3) = dN(1,a);  Bs(1, c+4) = N(a);
    Bs(2, c+3) = dN(2,a);  Bs(2, c+5) = N(a);
    Nm(:, c+1:c+5) = N(a) * eye(5);
end
end
