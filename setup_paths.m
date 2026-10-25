function root = setup_paths()
%SETUP_PATHS  Add all project folders to the MATLAB/Octave path.
%   Run once per session (every run_* script calls it automatically).
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'validation'));
addpath(fullfile(root, 'milestone1'));
addpath(fullfile(root, 'milestone2'));
addpath(fullfile(root, 'milestone3'));
addpath(fullfile(root, 'milestone4'));
addpath(fullfile(root, 'milestone5'));
if ~exist(fullfile(root, 'results', 'figures'), 'dir'), mkdir(fullfile(root, 'results', 'figures')); end
if ~exist(fullfile(root, 'results', 'data'), 'dir'),    mkdir(fullfile(root, 'results', 'data'));    end
if exist('OCTAVE_VERSION', 'builtin') && isempty(getenv('DISPLAY'))
    ws = warning('off', 'all'); graphics_toolkit('gnuplot'); warning(ws);  % headless Octave only
    warning('off', 'Octave:language-extension');
end
if exist('OCTAVE_VERSION', 'builtin')
    % Octave's default .mat format is a text format MATLAB cannot read;
    % always write MATLAB v7 binary files so results are portable
    save_default_options('-v7');
end
end
