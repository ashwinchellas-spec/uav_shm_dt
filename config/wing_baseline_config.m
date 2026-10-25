function cfg = wing_baseline_config()
%WING_BASELINE_CONFIG  Milestone 1 - problem definition of the composite UAV wing.
%
%   cfg = wing_baseline_config()  returns every modelling decision in one
%   struct, so later milestones (DOE, AI, digital twin) can perturb any
%   parameter programmatically.  All units SI.
%
%   Idealisation: the semi-wing is a cantilever sandwich plate clamped at
%   the fuselage (x = 0).  Carbon/epoxy skins on a structural PMI foam
%   core represent a foam-core composite UAV wing.  The core thickness is
%   chosen so the skins have the same chord-wise second moment as skins
%   following a NACA 0012 contour (see run_milestone1, section 2).

%% 1. Geometry (x: span, root -> tip;  y: chord, LE -> TE;  z: up)
cfg.name          = 'UAV-W1 baseline semi-wing';
cfg.geom.span     = 1.00;      % semi-span b              [m]
cfg.geom.c_root   = 0.25;      % root chord               [m]
cfg.geom.taper    = 1.0;       % c_tip / c_root           [-]
cfg.geom.sweep_le = 0.0;       % leading-edge sweep       [deg]
cfg.geom.airfoil  = 'NACA0012';
cfg.geom.tc       = 0.12;      % thickness-to-chord ratio [-]

%% 2. Aircraft / flight envelope
cfg.ac.mass   = 6.0;           % take-off mass             [kg]
cfg.ac.b_full = 2*cfg.geom.span;
cfg.ac.S      = 2*cfg.geom.span*cfg.geom.c_root*(1+cfg.geom.taper)/2;   % [m^2]
cfg.ac.cbar   = cfg.geom.c_root;
cfg.ac.rho    = 1.225;         % sea-level ISA density     [kg/m^3]
cfg.ac.CLmax  = 1.2;           % max lift coefficient      [-]
cfg.ac.CLmin  = -1.0;          % min (negative) lift coeff [-]
cfg.ac.Vc     = 25.0;          % design cruise speed       [m/s]
cfg.ac.Vd     = 32.0;          % design dive speed (>=1.25 Vc) [m/s]
cfg.ac.n_pos  = 3.8;           % manoeuvre limit (CS-VLA / CS-23 utility style)
cfg.ac.n_neg  = -1.5;
cfg.ac.Ude_c  = 15.24;         % derived gust at Vc (50 ft/s)  [m/s]
cfg.ac.Ude_d  = 7.62;          % derived gust at Vd (25 ft/s)  [m/s]
cfg.ac.SF     = 1.5;           % ultimate / limit safety factor
cfg.g         = 9.81;

%% 3. Materials and laminate (bottom -> top)
cfg.lam.face_mat    = 'T300_5208';
cfg.lam.face_angles = [0 45 -45 90];  % bottom skin, outer -> inner
cfg.lam.t_ply       = 0.125e-3;       % cured ply thickness  [m]
cfg.lam.core_mat    = 'Rohacell51WF';
cfg.lam.t_core      = 21.75e-3;       % core thickness       [m]
cfg.lam.kshear      = 'auto';         % energy-equivalent shear correction

%% 4. Finite-element model
cfg.fe.nx    = 40;             % elements along span (25 mm)
cfg.fe.ny    = 10;             % elements along chord (25 mm)
cfg.fe.shear = 'hr';           % Hellinger-Reissner mixed shear
cfg.fe.bc    = 'clamped_root';

%% 5. Analysis settings
cfg.an.n_modes   = 10;
cfg.an.zeta      = 0.02;       % modal damping ratio (Rayleigh, modes 1 & 5)
cfg.an.zeta_modes = [1 5];
cfg.an.tip_test_load = 50;     % ground-test tip load      [N]
cfg.an.gust_H_chords = 12.5;   % gust gradient H = 12.5 cbar (CS-25 style)
end
