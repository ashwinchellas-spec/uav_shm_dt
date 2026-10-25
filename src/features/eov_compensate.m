function [Xc, mdl] = eov_compensate(X, T, fitrows, mdl)
%EOV_COMPENSATE  Remove the temperature trend from every feature (regression).
%
%   [Xc, mdl] = eov_compensate(X, T, fitrows)   fit on rows 'fitrows'
%   Xc = eov_compensate(X, T, [], mdl)          apply a fitted model
%
%   For each feature j a linear model  x_j = a_j + b_j T  is fitted on
%   HEALTHY TRAINING data only (no leakage, no damage information), and the
%   residual  x_j - (a_j + b_j T)  is returned.  This is the classical
%   regression approach to environmental compensation (Sohn 2007; Peeters &
%   De Roeck 2001): what temperature explains is removed, what damage adds
%   remains.  Temperature itself is measured (noisy sensor).
if nargin < 4 || isempty(mdl)
    A = [ones(sum(fitrows),1), T(fitrows)];
    mdl.coef = A \ X(fitrows,:);
end
Xc = X - [ones(size(T)), T] * mdl.coef;
end
