function cd = discrepancy_cd2(X)
%DISCREPANCY_CD2  Centred L2-discrepancy of a design on [0,1]^k (Hickernell 1998).
%   Lower = more uniform.  Measures how far the empirical distribution of
%   the points deviates from the uniform distribution over all sub-boxes.
%   CD2^2 = (13/12)^k - (2/n) sum_i prod_j (1 + |z_ij|/2 - z_ij^2/2)
%           + (1/n^2) sum_i sum_l prod_j (1 + |z_ij|/2 + |z_lj|/2 - |x_ij - x_lj|/2),
%   z = x - 1/2.
[n, k] = size(X);
z = abs(X - 0.5);
t1 = (13/12)^k;
t2 = 2/n * sum(prod(1 + z/2 - z.^2/2, 2));
s3 = 0;
for i = 1:n
    s3 = s3 + sum(prod(1 + z(i,:)/2 + z/2 - abs(X(i,:) - X)/2, 2));
end
cd = sqrt(max(t1 - t2 + s3/n^2, 0));
end
