function m = simulate_measurement(db, i, sens, ref, noise)
%SIMULATE_MEASUREMENT  What the real sensor system would report for sample i.
%
%   m = simulate_measurement(db, i, sens, ref, noise)
%   Turns the noise-free simulation of database row i into a realistic
%   measurement record.  Error sources (noise struct, all optional; set a
%   field to 0 for a noise-free measurement):
%     .T      temperature sensor, std [C]                      (0.5)
%     .f      modal-identification error of frequencies, rel.  (5e-4)
%     .phi    mode-shape identification noise, x max|phi|      (0.01)
%     .ms     dynamic FBG modal-strain noise, x max|eps|       (0.02)
%     .eps    static FBG strain noise [strain]                 (1e-6)
%     .gain   FBG gauge-factor (gain) error per grating, rel.  (0.005)
%     .load   tip-test load uncertainty, rel.                  (0.01)
%     .field  distributed-fibre (dense map) strain noise       (2e-6)
%
%   Modal data: z-accelerometers only observe flexural modes, so the
%   in-plane mode is not identified.  Identified modes are PAIRED with the
%   reference modes by maximum MAC (as in real operational modal analysis),
%   never by their index, and their sign is aligned with the reference.
%
%   m.T, m.f (1x5), m.phi (n_acc x 5), m.ms (n_fbg x 5), m.pair_mac (1x5)
%   m.eps  (n_fbg x 2)  static FBG strain, bottom and TOP skin (both faces
%                       are instrumented: a sensor on one skin barely sees
%                       damage in the other skin of a sandwich)
%   m.field (ne x 3 x 2) dense static strain [ex ey gxy] bottom / top
%                       (full-field measurement, e.g. DIC in a ground test)
%   m.mfield (ne x 3 x 5) dense bottom-skin modal strains of the 5 modes
def = struct('T',0.5,'f',5e-4,'phi',0.01,'ms',0.02,'eps',1e-6,'gain',0.005,'load',0.01,'field',2e-6);
if nargin < 5 || isempty(noise), noise = def; end
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(noise, fn{k}), noise.(fn{k}) = def.(fn{k}); end, end

ia = sens.acc_node; ie = sens.fbg_elem;
m.T = db.P(i, strcmp(db.var_names, 'T')) + noise.T*randn;

cand = find(db.mode_type(i,:) ~= 4);               % observable (flexural) modes
fc  = db.freq(i, cand) .* (1 + noise.f*randn(1, numel(cand)));
Pc  = double(db.Wmodes(ia, cand, i));
Pc  = Pc + noise.phi*max(abs(Pc), [], 1).*randn(size(Pc));
Mc  = double(squeeze(db.Emodes(ie, 1, cand, i)));
Mc  = Mc + noise.ms*max(abs(Mc), [], 1).*randn(size(Mc));
Pr  = double(ref.o.Wmodes(ia, ref.flex));
nr = numel(ref.flex);
mac = ((Pr'*Pc).^2) ./ (sum(Pr.^2,1)' * sum(Pc.^2,1));
m.f = nan(1,nr); m.phi = nan(numel(ia), nr); m.ms = nan(numel(ie), nr); m.pair_mac = nan(1,nr);
used = false(1, numel(cand)); kk = zeros(1, nr);
for j = 1:nr
    mj = mac(j,:); mj(used) = -inf;
    [m.pair_mac(j), k] = max(mj); used(k) = true; kk(j) = k;
    s = sign(Pr(:,j)'*Pc(:,k)); if s == 0, s = 1; end
    m.f(j) = fc(k); m.phi(:,j) = s*Pc(:,k); m.ms(:,j) = s*Mc(:,k);
end

% static tip test (L2): gain error per grating, load uncertainty, noise
lf = 1 + noise.load*randn;
e = [double(db.eps_bot(ie, 1, 2, i)), double(db.eps_top(ie, 1, 2, i))];
m.eps = e .* (1 + noise.gain*randn(size(e))) * lf + noise.eps*randn(size(e));
% dense full-field static strain maps (both skins, 3 components)
ef = cat(3, double(db.eps_bot(:, :, 2, i)), double(db.eps_top(:, :, 2, i)));
m.field = ef * lf + noise.field*randn(size(ef));
% dense modal strain maps (bottom skin, 3 components) of the paired modes
Ef = double(db.Emodes(:, :, cand, i));
m.mfield = nan(size(Ef,1), 3, nr);
for j = 1:nr
    k = kk(j);
    s = sign(Pr(:,j)'*Pc(:,k)); if s == 0, s = 1; end
    E = Ef(:,:,k);
    m.mfield(:,:,j) = s*(E + noise.ms*max(abs(E(:)))*randn(size(E)));
end
end
