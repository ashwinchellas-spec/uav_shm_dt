function auc = roc_auc(score, label)
%ROC_AUC  Area under the ROC curve (Mann-Whitney U form, ties averaged).
%   auc = roc_auc(score, label)   label: 1 = positive (damaged), 0 = negative
%   AUC = P(score of a random positive > score of a random negative).
%   0.5 = chance, 1 = perfect separation.
score = score(:); label = logical(label(:));
n1 = sum(label); n0 = sum(~label);
[~, o] = sort(score); r = zeros(size(score));
r(o) = 1:numel(score);
% average ranks of ties
[us, ~, g] = unique(score);
if numel(us) < numel(score)
    for k = 1:numel(us)
        q = g == k; if sum(q) > 1, r(q) = mean(r(q)); end
    end
end
auc = (sum(r(label)) - n1*(n1+1)/2) / (n1*n0);
end
