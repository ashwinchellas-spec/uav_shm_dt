function B = padarray_rep(A)
%PADARRAY_REP  Pad a matrix by one replicated row/column on every side
%   (no Image Processing Toolbox needed); used before 3x3 smoothing.
B = [A(1,1) A(1,:) A(1,end); A(:,1) A A(:,end); A(end,1) A(end,:) A(end,end)];
end
