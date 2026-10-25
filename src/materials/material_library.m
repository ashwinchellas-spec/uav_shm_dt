function mat = material_library(name)
%MATERIAL_LIBRARY  Lamina / core property database (SI units).
%
%   mat = material_library(name) returns a struct with the orthotropic
%   lamina properties used by the FSDT solver:
%       E1, E2      Young's moduli along / across fibre          [Pa]
%       G12         in-plane shear modulus                       [Pa]
%       G13, G23    transverse shear moduli                      [Pa]
%       nu12        major Poisson ratio                          [-]
%       rho         density                                      [kg/m^3]
%       Xt,Xc,Yt,Yc,S12,S23  strengths (positive numbers)        [Pa]
%       type        'lamina' | 'core' | 'isotropic'
%
%   Available names:
%       'T300_5208'   carbon/epoxy (Tsai, 1980)  - baseline wing skins
%       'AS4_3501_6'  carbon/epoxy (Daniel & Ishai) - DOE alternative
%       'Eglass_epoxy' E-glass/epoxy (Scotchply 1002) - DOE alternative
%       'Rohacell51WF' PMI structural foam core - baseline wing core
%       'Al7075'      aluminium (isotropic, validation)
%       'Reddy_E25'   E1/E2=25 benchmark material (Reddy, Pagano)
%       'Reddy_E40'   E1/E2=40 benchmark material (Reddy, Noor)
%       'iso_unit'    isotropic E=1, nu=0.3 (Leissa / locking tests)
%
%   Transverse shear modulus G23 of a transversely isotropic lamina is
%   G23 = E2 / (2 (1 + nu23)).  nu23 ~ 0.49 is used for carbon/epoxy.

switch lower(name)
    case 't300_5208'
        mat = lamina(181e9, 10.3e9, 7.17e9, 7.17e9, 10.3e9/(2*1.49), 0.28, 1600, ...
                     1500e6, 1500e6, 40e6, 246e6, 68e6, 40e6);
    case 'as4_3501_6'
        mat = lamina(147e9, 10.3e9, 7.0e9, 7.0e9, 3.7e9, 0.27, 1580, ...
                     2280e6, 1725e6, 57e6, 228e6, 76e6, 50e6);
    case 'eglass_epoxy'
        mat = lamina(38.6e9, 8.27e9, 4.14e9, 4.14e9, 8.27e9/(2*1.4), 0.26, 1800, ...
                     1062e6, 610e6, 31e6, 118e6, 72e6, 40e6);
    case 'rohacell51wf'
        % PMI foam: E = 75 MPa, G = 24 MPa (datasheet). A foam is not a
        % perfect isotropic solid (E/2G-1 > 0.5), so it is stored as a
        % quasi-isotropic orthotropic layer with independent G.
        mat = lamina(75e6, 75e6, 24e6, 24e6, 24e6, 0.30, 52, ...
                     1.6e6, 0.8e6, 1.6e6, 0.8e6, 0.8e6, 0.8e6);
        mat.type = 'core';
    case 'al7075'
        E = 71.7e9; nu = 0.33; G = E/(2*(1+nu));
        mat = lamina(E, E, G, G, G, nu, 2810, 503e6, 503e6, 503e6, 503e6, 331e6, 331e6);
        mat.type = 'isotropic';
    case 'reddy_e25'
        E2 = 1e9;   % only ratios matter for normalised benchmarks
        mat = lamina(25*E2, E2, 0.5*E2, 0.5*E2, 0.2*E2, 0.25, 1, 1,1,1,1,1,1);
    case 'reddy_e40'
        E2 = 1e9;
        mat = lamina(40*E2, E2, 0.6*E2, 0.6*E2, 0.5*E2, 0.25, 1, 1,1,1,1,1,1);
    case 'iso_unit'
        E = 1; nu = 0.3; G = E/(2*(1+nu));
        mat = lamina(E, E, G, G, G, nu, 1, 1,1,1,1,1,1);
        mat.type = 'isotropic';
    otherwise
        error('material_library: unknown material "%s"', name);
end
mat.name = name;
end

function m = lamina(E1,E2,G12,G13,G23,nu12,rho,Xt,Xc,Yt,Yc,S12,S23)
m.E1 = E1; m.E2 = E2; m.G12 = G12; m.G13 = G13; m.G23 = G23;
m.nu12 = nu12; m.rho = rho;
m.Xt = Xt; m.Xc = Xc; m.Yt = Yt; m.Yc = Yc; m.S12 = S12; m.S23 = S23;
m.type = 'lamina';
end
