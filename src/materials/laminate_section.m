function sec = laminate_section(layup, kopt)
%LAMINATE_SECTION  FSDT section properties of a laminate (CLT + shear).
%
%   sec = laminate_section(layup)
%   sec = laminate_section(layup, kopt)
%
%   layup.mat   : 1xN struct array of materials (bottom -> top)
%   layup.theta : 1xN ply angles [deg]
%   layup.t     : 1xN ply thicknesses [m]
%   kopt        : shear correction - scalar (e.g. 5/6, applied to both
%                 directions) or 'auto' (energy-equivalent factors computed
%                 from the actual stacking, see shear_correction.m).
%                 Default 5/6.
%
%   Returned fields
%     A,B,D (3x3)   membrane, coupling, bending stiffness
%     As    (2x2)   transverse shear stiffness  [Q_x Q_y]' = As [g_xz g_yz]'
%     I0,I1,I2      mass moments  (int rho z^k dz)
%     h, z          total thickness and ply interfaces (z(1) = -h/2)
%     Qb{k},Qsb{k},Te{k},Qm{k}  transformed / material ply stiffness
%     kx, ky        shear correction factors used
%     tprof         through-thickness shear-stress profile (stress recovery)

if nargin < 2 || isempty(kopt), kopt = 5/6; end

N = numel(layup.theta);
t = layup.t(:)';
h = sum(t);
z = [-h/2, -h/2 + cumsum(t)];

A = zeros(3); B = zeros(3); D = zeros(3); Asr = zeros(2);
I0 = 0; I1 = 0; I2 = 0;
Qb = cell(1,N); Qsb = cell(1,N); Te = cell(1,N); Qm = cell(1,N);
for k = 1:N
    [Q, Qs] = ply_stiffness(layup.mat(k));
    Qm{k} = Q;
    [Qb{k}, Qsb{k}, Te{k}] = transform_ply(Q, Qs, layup.theta(k));
    z1 = z(k); z2 = z(k+1);
    A   = A   + Qb{k} * (z2 - z1);
    B   = B   + Qb{k} * (z2^2 - z1^2) / 2;
    D   = D   + Qb{k} * (z2^3 - z1^3) / 3;
    Asr = Asr + Qsb{k} * (z2 - z1);
    r = layup.mat(k).rho;
    I0 = I0 + r * (z2 - z1);
    I1 = I1 + r * (z2^2 - z1^2) / 2;
    I2 = I2 + r * (z2^3 - z1^3) / 3;
end

sec.layup = layup;
sec.h = h; sec.z = z; sec.nply = N;
sec.A = A; sec.B = B; sec.D = D;
sec.Qb = Qb; sec.Qsb = Qsb; sec.Te = Te; sec.Qm = Qm;
sec.I0 = I0; sec.I1 = I1; sec.I2 = I2;

% Shear correction factors and through-thickness shear-stress profiles
sc = shear_correction(sec);
if ischar(kopt) && strcmpi(kopt, 'auto')
    kx = sc.kx; ky = sc.ky;
else
    kx = kopt; ky = kopt;
end
sec.kx = kx; sec.ky = ky;
sec.k_energy = [sc.kx sc.ky];
sec.As = [kx*Asr(1,1), sqrt(kx*ky)*Asr(1,2); sqrt(kx*ky)*Asr(2,1), ky*Asr(2,2)];
sec.As_raw = Asr;
sec.tprof = sc.tprof;
end
