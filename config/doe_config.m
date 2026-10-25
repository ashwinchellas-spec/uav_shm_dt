function spec = doe_config(tier)
%DOE_CONFIG  Milestone 3 - Design-of-Experiments variable spaces.
%
%   spec = doe_config('shm')     main SHM database (fixed wing design)
%   spec = doe_config('design')  design-space database (healthy, varied design)
%
%   Each continuous variable: name, lower, upper, scale ('lin'|'log'),
%   unit, group.  Values are drawn by maximin Latin Hypercube on [0,1]^k
%   and mapped  x = lo + u (hi - lo);  for 'log' variables lo and hi are
%   log10 EXPONENTS and x = 10^(lo + u (hi - lo)).
%
%   TIER 'shm' - the digital twin of ONE wing in service
%     * manufacturing scatter of the nominal design (material, thickness,
%       ply alignment) - the "aleatory" variability every real wing has
%     * environmental & operational variability (temperature, root fixity)
%       - the confounders that make real SHM hard
%     * damage scenarios: 6 damage classes + healthy, stratified (equal
%       counts), each with location, size, severity, face
%     Load magnitude is NOT a DOE variable: the model is linear, so the
%     response to any load factor n is n x (stored unit-load response).
%     Unit 1 g, tip-load and gust responses are stored instead.
%
%   TIER 'design' - generalisation over designs (healthy wings)
%     ply orientations, ply / core thickness (thickness ratio), face
%     material (categorical), core grade, planform taper & sweep, root
%     fixity - used for the parametric surrogate / PINN in Milestone 5.

spec.tier = tier;
switch lower(tier)
    case 'shm'
        spec.N_per_class = 400;
        spec.classes = {'healthy', 'debond', 'skin_matrix', 'skin_fiber', ...
                        'core_crush', 'delam', 'root_loosening'};
        spec.N = spec.N_per_class * numel(spec.classes);
        spec.seed = 2024;
        v = {
        % name      lo      hi     scale  unit     group
        'E1'        0.95    1.05   'lin'  '-'      'material'
        'E2'        0.92    1.08   'lin'  '-'      'material'
        'G12'       0.92    1.08   'lin'  '-'      'material'
        'Gc'        0.85    1.15   'lin'  '-'      'material'
        'Ec'        0.85    1.15   'lin'  '-'      'material'
        'rho_f'     0.97    1.03   'lin'  '-'      'material'
        'rho_c'     0.90    1.10   'lin'  '-'      'material'
        'tply'      0.95    1.05   'lin'  '-'      'geometry'
        'tcore'     0.97    1.03   'lin'  '-'      'geometry'
        'dtheta'   -2.0     2.0    'lin'  'deg'    'geometry'
        'T'        -10      40     'lin'  'degC'   'environment'
        'kroot'     2.3     3.7    'log'  '-'      'bc'          % 10^[...] = 200 ... 5000
        'dmg_x'     0.03    0.97   'lin'  'x/L'    'damage'
        'dmg_y'     0.12    0.88   'lin'  'y/c'    'damage'
        'dmg_r'     0.015   0.060  'lin'  'm'      'damage'
        'dmg_ar'    0.6     1.4    'lin'  '-'      'damage'
        'dmg_D'     0.10    0.90   'lin'  '-'      'damage'
        'dmg_face'  0       1      'lin'  '-'      'damage'
        };
    case 'design'
        spec.N_per_class = 400;
        spec.classes = {'healthy'};
        spec.N = 400;
        spec.seed = 4048;
        v = {
        'th1'      -90      90     'lin'  'deg'    'layup'
        'th2'      -90      90     'lin'  'deg'    'layup'
        'th3'      -90      90     'lin'  'deg'    'layup'
        'th4'      -90      90     'lin'  'deg'    'layup'
        'tply_mm'   0.10    0.20   'lin'  'mm'     'geometry'
        'tcore_mm'  8       30     'lin'  'mm'     'geometry'
        'core_grade' -0.30103 0.30103 'log' '-'     'material'   % 10^[...] = 0.5 ... 2
        'face_mat'  0       1      'lin'  '-'      'material'
        'taper'     0.5     1.0    'lin'  '-'      'geometry'
        'sweep'     0       20     'lin'  'deg'    'geometry'
        'kroot'     2.3     3.7    'log'  '-'      'bc'
        };
        spec.face_mats = {'T300_5208', 'AS4_3501_6', 'Eglass_epoxy'};
    otherwise
        error('doe_config: unknown tier');
end
spec.vars = cell2struct(v, {'name','lo','hi','scale','unit','group'}, 2);
spec.k = numel(spec.vars);
spec.n_modes = 8;
spec.gust_tg = [];                 % [] -> design gust from the flight envelope
spec.zeta = 0.02;                  % modal damping at modes 1 and 5 (Rayleigh)
spec.x_off = 0.05;                 % strength read-out station (root singularity)
end
