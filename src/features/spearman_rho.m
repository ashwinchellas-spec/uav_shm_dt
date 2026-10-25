function rho = spearman_rho(a, b)
%SPEARMAN_RHO  Rank correlation (monotonic association), no toolbox needed.
[~, i] = sort(a(:)); ra = zeros(size(i)); ra(i) = 1:numel(i);
[~, j] = sort(b(:)); rb = zeros(size(j)); rb(j) = 1:numel(j);
c = corrcoef(ra, rb); rho = c(1,2);
end
