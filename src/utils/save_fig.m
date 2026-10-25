function save_fig(fig, name)
%SAVE_FIG  Save a figure as PNG in results/figures (MATLAB & Octave).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
f = fullfile(root, 'results', 'figures', [name '.png']);
set(fig, 'PaperPositionMode', 'auto');
print(fig, '-dpng', '-r130', f);
close(fig);
end
