function split = split_by_sample(cls, fr, seed)
%SPLIT_BY_SAMPLE  Stratified train / validation / test split of SAMPLES.
%   split = split_by_sample(cls, fr, seed)    fr = [0.70 0.15 0.15]
%   split(i) = 1 train, 2 validation, 3 test, for every simulated wing i.
%   All noisy realisations of one wing must end up in the SAME subset;
%   otherwise the model sees (almost) the test wing during training and
%   the test score is optimistic (data leakage).  Stratified so every
%   class has the same proportions in every subset.
rand('state', seed); %#ok<RAND>
split = zeros(numel(cls), 1);
for c = unique(cls(:))'
    ii = find(cls == c);
    ii = ii(randperm(numel(ii)));
    n1 = round(fr(1)*numel(ii)); n2 = round(fr(2)*numel(ii));
    split(ii(1:n1)) = 1; split(ii(n1+1:n1+n2)) = 2; split(ii(n1+n2+1:end)) = 3;
end
end
