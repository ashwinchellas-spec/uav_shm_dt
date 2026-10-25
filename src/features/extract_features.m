function [x, names, groups, maps] = extract_features(m, r)
%EXTRACT_FEATURES  Physics-informed damage features from one measurement.
%
%   [x, names, groups, maps] = extract_features(m, r)
%   m : measurement (simulate_measurement);  r : noise-free REFERENCE
%       measurement of the baseline (healthy, as-designed) wing.
%
%   Groups (each is designed against a known physical effect):
%   'freq'    f_k [Hz] of the 5 flexural modes (MAC-paired)
%   'ratio'   f_k / f_1, k = 2..5.  A uniform change of stiffness or mass
%             (temperature, batch) scales ALL frequencies alike and cancels
%             in the ratio; local damage changes them unequally.
%   'mac'     1 - MAC(phi, phi_ref) per mode (mode-shape change)
%   'uls'     uniform-load-surface residual at the accelerometers.  Modal
%             flexibility G = Phi Omega^-2 Phi' (mass-normalised modes) is
%             dominated by the low modes and is a direct stiffness measure;
%             u = G 1.  After removing the best global scale a,
%                 r = (u/a - u_ref) / max|u_ref|          (Pandey & Biswas 1994)
%   'strain'  static tip-test strain SHAPE residual at the FBGs on BOTH
%             skins (60 values), same global-scale removal (Milestone 3:
%             EOV acts almost as a global scaling of the strain field)
%   'msedi'   log of the modal-strain-energy damage index per FBG
%             (Stubbs / Cornwell, sensor form):
%                 F_sk = eps_sk^2 / sum_s eps_sk^2      (energy fraction)
%                 beta_s = sum_k F_sk / sum_k F_sk^ref
%             fractions are scale-free -> insensitive to global EOV.
%             A noise floor lambda = 0.5 mean(F_ref) is added to numerator
%             and denominator: where the reference modal strain is ~0
%             (free tip, nodal lines) noise would otherwise explode the
%             ratio (the role of the '+1' in Stubbs' original index).
%   'temp'    measured temperature (for EOV compensation by the models)
%
%   maps (dense full-field data, for localisation and CNN models):
%     maps.strain : static strain-shape residual magnitude per element,
%                   both skins, 3 components (ne x 1)
%     maps.msedi  : log MSE damage index per element from the 3-component
%                   bottom-skin modal strains (ne x 1)
%   Before localisation both maps are z-scored against the HEALTHY
%   population (per element), see run_milestone4.
sc = @(u, u0) (u0'*u) / (u0'*u0);                 % best global scale factor

f = m.f(:)';
x_freq = f;
x_ratio = f(2:end) / f(1);
mac = (sum(m.phi.*r.phi,1).^2) ./ (sum(m.phi.^2,1).*sum(r.phi.^2,1));
x_mac = 1 - mac;

w = 2*pi*f; w0 = 2*pi*r.f;
u  = m.phi * ((m.phi' * ones(size(m.phi,1),1)) ./ w(:).^2);
u0 = r.phi * ((r.phi' * ones(size(r.phi,1),1)) ./ w0(:).^2);
x_uls = ((u / sc(u, u0)) - u0)' / max(abs(u0));

e = m.eps(:); e0 = r.eps(:);
x_strain = ((e / sc(e, e0)) - e0)' / max(abs(e0));

F  = sum(m.ms.^2 ./ sum(m.ms.^2, 1), 2);
F0 = sum(r.ms.^2 ./ sum(r.ms.^2, 1), 2);
lam = 0.5*mean(F0);
x_msedi = log((F + lam) ./ (F0 + lam))';

x_temp = m.T;

x = [x_freq, x_ratio, x_mac, x_uls, x_strain, x_msedi, x_temp];
na = numel(x_uls); ns = numel(x_strain);
names = [arrayfun(@(k) sprintf('f%d', k), 1:5, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('f%d/f1', k), 2:5, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('1-MAC%d', k), 1:5, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('ULS_a%02d', k), 1:na, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('epsB_s%02d', k), 1:ns/2, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('epsT_s%02d', k), 1:ns/2, 'UniformOutput', false), ...
         arrayfun(@(k) sprintf('MSEDI_s%02d', k), 1:numel(x_msedi), 'UniformOutput', false), {'T'}];
groups = [repmat({'freq'},1,5), repmat({'ratio'},1,4), repmat({'mac'},1,5), ...
          repmat({'uls'},1,na), repmat({'strain'},1,ns), repmat({'msedi'},1,numel(x_msedi)), {'temp'}];

if nargout > 3
    ef = m.field(:); ef0 = r.field(:);
    d = reshape(((ef / sc(ef, ef0)) - ef0) / max(abs(ef0)), [], 6);
    maps.strain = sqrt(sum(d.^2, 2));
    Um  = squeeze(sum(m.mfield.^2, 2)); Um0 = squeeze(sum(r.mfield.^2, 2));
    Fm  = sum(Um ./ sum(Um, 1), 2); Fm0 = sum(Um0 ./ sum(Um0, 1), 2);
    lm = 0.5*mean(Fm0);
    maps.msedi = log((Fm + lm) ./ (Fm0 + lm));
end
end
