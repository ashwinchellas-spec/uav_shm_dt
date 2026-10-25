function [p, w] = gauss_1d(n)
%GAUSS_1D  Gauss-Legendre points and weights on [-1,1] (n = 1..6).
switch n
    case 1, p = 0; w = 2;
    case 2, p = [-1 1]/sqrt(3); w = [1 1];
    case 3, p = [-sqrt(3/5) 0 sqrt(3/5)]; w = [5 8 5]/9;
    otherwise
        % Golub-Welsch
        k = 1:n-1; bet = k ./ sqrt(4*k.^2 - 1);
        [V, L] = eig(diag(bet,1) + diag(bet,-1));
        [p, i] = sort(diag(L)'); w = 2 * V(1,i).^2;
end
end
