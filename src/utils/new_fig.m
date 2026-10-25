function fig = new_fig(w, h)
%NEW_FIG  Invisible figure of w x h pixels (for batch / headless runs).
if nargin < 1, w = 900; end
if nargin < 2, h = 600; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 w h]);
end
