function [Xin, info] = novelty_inputs(X, groups, href, info)
%NOVELTY_INPUTS  Healthy-referenced, whitened network inputs.
%
%   [Xin, info] = novelty_inputs(X, groups, href)    fit on rows href
%   Xin = novelty_inputs(X, groups, [], info)        apply a stored fit
%
%   Damage is a DEVIATION from the healthy population, and that population
%   is strongly correlated by the operational/environmental variability
%   (temperature, root fixity, mass, material scatter move many features
%   together).  The inputs given to the network are therefore
%     1) the Mahalanobis-whitened residual  z = S_h^(-1/2) (x - mu_h)
%        (healthy training rows only; S_h shrunk 5 % towards its diagonal),
%        so that the EOV directions are already normalised away, and
%     2) log(1 + D_g^2/n_g) for every feature group g (frequency, ratio,
%        MAC, ULS, strain, MSE index): the classical novelty indices of
%        Milestone 4, one per physical family.
%   The network then only has to LEARN what statistics cannot do: which
%   pattern of deviations belongs to which damage type, where and how big.
if nargin < 4 || isempty(info)
    Xh = X(href, :);
    info.mu = mean(Xh, 1);
    S = cov(Xh); a = 0.05; S = (1 - a)*S + a*diag(diag(S));
    [V, E] = eig((S + S')/2); e = max(diag(E), 1e-10*max(diag(E)));
    info.W = V * diag(1 ./ sqrt(e)) * V';           % symmetric S^(-1/2)
    info.gnames = unique(groups, 'stable');
    info.gnames = info.gnames(~strcmp(info.gnames, 'temp'));   % temperature is EOV, not damage
    for k = 1:numel(info.gnames)
        c = strcmp(groups, info.gnames{k});
        [~, info.mdl{k}] = novelty_mahal(Xh(:, c), Xh(1, c));
    end
end
Zw = (X - info.mu) * info.W;
G = zeros(size(X,1), numel(info.gnames));
for k = 1:numel(info.gnames)
    c = strcmp(groups, info.gnames{k});
    G(:,k) = log(1 + novelty_mahal([], X(:, c), info.mdl{k}) / sum(c));
end
if ~isfield(info, 'gmu')
    info.gmu = mean(G(href,:), 1); info.gsd = std(G(href,:), 0, 1);
end
Xin = [Zw, (G - info.gmu) ./ info.gsd];
end
