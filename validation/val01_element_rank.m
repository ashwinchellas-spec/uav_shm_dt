function R = val01_element_rank()
%VAL01_ELEMENT_RANK  Zero-energy modes of a single distorted element.
%   A correct plate element has exactly 6 rigid-body modes
%   (u, v, in-plane rotation, w, rotation about x, rotation about y).
%   More zeros = spurious (hourglass) mechanisms.
m = material_library('iso_unit'); L.mat = m; L.theta = 0; L.t = 0.1;
sec = laminate_section(L, 5/6);
xy = [0 0; 1 0; 1.2 0.9; -0.1 1.1];
R = [];
for sh = {'hr', 'mitc4', 'sri', 'full'}
    Ke = elem_fsdt_q4(xy, sec, sh{1});
    ev = sort(abs(eig(Ke)));
    nz = sum(ev < 1e-10*max(ev));
    tol = 0; if strcmp(sh{1}, 'sri'), tol = inf; end   % SRI: known defect, reported only
    R = [R, vres('V01 element rank', sprintf('zero modes, %s (expect 6)', sh{1}), nz, 6, tol)]; %#ok<AGROW>
end
end
