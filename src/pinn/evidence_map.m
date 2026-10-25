function E = evidence_map(F, href, mode, xy, ell)
%EVIDENCE_MAP  Physics-based damage-evidence map used as a prior by PI-Net.
%
%   E = evidence_map(F, href, mode, xy, ell)
%   F    : Milestone 4 feature set;  href : healthy training rows
%   mode : 'fbg'   - 60 FBG sensors only (realistic network): the
%                    healthy-z-scored MSE index (bottom FBGs) and |strain
%                    residual| (bottom + top FBGs) of Milestone 4 are summed per
%                    grating position and spread to the elements with a
%                    Gaussian kernel of length ell (default 60 mm);
%          'dense' - full-field map (Milestone 4 dense z-map: modal + static
%                    two-skin evidence on every element), e.g. from DIC or a
%                    distributed (Rayleigh-backscatter) fibre;
%          'none'  - zeros (no prior).
%   xy   : element centroids (ne x 2).
%   E (n x ne), scaled to unit standard deviation over healthy rows so the
%   learnable gain of the map head starts at a sensible value.
%
%   The map head of PI-Net then computes  q = softmax(W h + b + a E): the
%   network REFINES a physics-based evidence map instead of learning the
%   geometry of the wing from scratch (residual / prior learning).
if nargin < 5, ell = 0.06; end
n = size(F.X, 1); ne = size(xy, 1);
switch mode
    case 'none'
        E = zeros(n, ne); return
    case 'dense'
        E = double(F.map_z);
    case 'fbg'
        Xc = double(F.X_comp); g = F.groups;
        im = find(strcmp(g, 'msedi')); is = find(strcmp(g, 'strain'));
        zz = @(A) (A - mean(A(href,:), 1)) ./ std(A(href,:), 0, 1);
        zs = zz(Xc(:, im)) + zz(abs(Xc(:, is(1:30)))) + zz(abs(Xc(:, is(31:60))));
        pos = F.sensors.fbg_pos;
        d2 = (xy(:,1) - pos(:,1)').^2 + (xy(:,2) - pos(:,2)').^2;
        K = exp(-d2 / (2*ell^2)); K = K ./ sum(K, 2);
        E = zs * K';
    otherwise
        error('evidence_map: unknown mode');
end
E = (E - mean(E(href,:), 1)) / std(reshape(E(href,:), [], 1));
end
