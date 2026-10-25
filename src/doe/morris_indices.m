function S = morris_indices(U, Y, delta, k)
%MORRIS_INDICES  Elementary effects -> mu* (importance) and sigma (interaction).
%   S = morris_indices(U, Y, delta, k)
%   U : r(k+1) x k design from morris_design;  Y : r(k+1) x m outputs
%   EE_i = [y(x + delta e_i) - y(x)] / delta     (in normalised units)
%   mu*_i  = mean |EE_i|   (Campolongo et al. 2007) - overall influence
%   sigma_i = std(EE_i)    - non-linearity / interaction with other factors
m = size(Y, 2); r = size(U,1)/(k+1);
EE = zeros(r, k, m);
for t = 1:r
    idx = (t-1)*(k+1) + (1:k+1);
    Ut = U(idx,:); Yt = Y(idx,:);
    for s = 1:k
        dU = Ut(s+1,:) - Ut(s,:);
        i = find(abs(dU) > 1e-12, 1);
        EE(t, i, :) = reshape((Yt(s+1,:) - Yt(s,:)) / dU(i), 1, 1, m);
    end
end
S.mustar = reshape(mean(abs(EE), 1), k, m);
S.mu     = reshape(mean(EE, 1), k, m);
S.sigma  = reshape(std(EE, 0, 1), k, m);
S.EE = EE;
end
