function xy = localize_damage(idx, pos, k)
%LOCALIZE_DAMAGE  Damage location estimate from a damage-index map.
%   xy = localize_damage(idx, pos, k)
%   idx : n x 1 damage index (larger = more likely damaged)
%   pos : n x 2 coordinates of the index points (elements or sensors)
%   k   : number of strongest points used (default 3)
%   xy  : index-weighted centroid of the k largest positive values.
if nargin < 3, k = 3; end
v = idx(:); v(~isfinite(v)) = -inf;
[s, o] = sort(v, 'descend');
o = o(1:k); s = max(s(1:k), 0);
if sum(s) <= 0, xy = pos(o(1),:); return; end
xy = (s' * pos(o,:)) / sum(s);
end
