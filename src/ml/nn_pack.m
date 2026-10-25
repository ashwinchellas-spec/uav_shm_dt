function [theta, shapes] = nn_pack(P)
%NN_PACK  Flatten a cell array of parameter matrices into one column vector.
%   [theta, shapes] = nn_pack(P);   P = nn_unpack(theta, shapes)
%   One vector makes the optimiser (Adam) and the gradient check generic.
shapes = cellfun(@size, P, 'UniformOutput', false);
theta = cell2mat(cellfun(@(A) A(:), P(:), 'UniformOutput', false));
end
