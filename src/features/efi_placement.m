function [sel, Ed_hist] = efi_placement(Phi, cand, ns)
%EFI_PLACEMENT  Effective Independence sensor placement (Kammer 1991).
%
%   [sel, Ed_hist] = efi_placement(Phi, cand, ns)
%   Phi  : n_cand x m target mode shapes at the candidate locations
%   cand : candidate ids (returned subset is taken from these)
%   ns   : number of sensors to keep (ns >= m)
%
%   The Fisher information matrix of the modal coordinates is
%   Q = Phi' Phi.  The contribution of candidate i to its rank / volume is
%       E_D,i = [Phi (Phi' Phi)^-1 Phi']_ii       (0 <= E_D <= 1)
%   Candidates are removed one at a time, always the one with the smallest
%   E_D, until ns remain.  The kept set maximises det(Q) approximately, so
%   the target modes stay linearly independent and observable in noise -
%   exactly what is needed to identify mode shapes (MAC, flexibility).
cand = cand(:); idx = (1:numel(cand))';
Ed_hist = [];
while numel(idx) > ns
    P = Phi(idx,:);
    Ed = sum((P / (P'*P)) .* P, 2);
    [mn, k] = min(Ed);
    Ed_hist(end+1) = mn; %#ok<AGROW>
    idx(k) = [];
end
sel = cand(idx);
end
