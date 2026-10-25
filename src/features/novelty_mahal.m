function [score, mdl] = novelty_mahal(Xtrain_healthy, X, mdl)
%NOVELTY_MAHAL  Mahalanobis-distance novelty detector (outlier analysis).
%
%   [score, mdl] = novelty_mahal(Xh, X)      fit on healthy data Xh, score X
%   score = novelty_mahal([], X, mdl)        score with an existing model
%
%   D^2(x) = (x - mu)' S^-1 (x - mu),  mu, S from HEALTHY training data only
%   (Worden, Manson & Fieller 2000).  The covariance is shrunk towards its
%   diagonal (Ledoit-Wolf style, fixed weight 0.05) to stay well conditioned
%   when features are many and correlated.  A sample is flagged as damaged
%   when D^2 exceeds a threshold set on healthy data (e.g. the 99th
%   percentile -> 1 % false alarms).
if nargin < 3 || isempty(mdl)
    mdl.mu = mean(Xtrain_healthy, 1);
    S = cov(Xtrain_healthy);
    a = 0.05;
    S = (1 - a)*S + a*diag(diag(S));
    mdl.Sinv = inv(S + 1e-12*trace(S)/size(S,1)*eye(size(S)));
end
Z = X - mdl.mu;
score = sum((Z * mdl.Sinv) .* Z, 2);
end
