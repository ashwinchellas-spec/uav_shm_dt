function P = nn_unpack(theta, shapes)
%NN_UNPACK  Inverse of nn_pack.
P = cell(size(shapes)); k = 0;
for i = 1:numel(shapes)
    n = prod(shapes{i});
    P{i} = reshape(theta(k+1:k+n), shapes{i});
    k = k + n;
end
end
