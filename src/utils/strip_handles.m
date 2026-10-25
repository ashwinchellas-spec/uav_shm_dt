function s = strip_handles(s)
%STRIP_HANDLES  Remove function-handle fields (recursively) before saving.
%   MATLAB v7 .mat files written by Octave cannot hold function handles,
%   and handles are code rather than results - they can always be rebuilt
%   (e.g. env = flight_envelope(cfg.ac)).  Struct arrays are handled too.
if isstruct(s)
    fn = fieldnames(s);
    for k = 1:numel(fn)
        if isa(s(1).(fn{k}), 'function_handle')
            s = rmfield(s, fn{k});
        else
            for j = 1:numel(s), s(j).(fn{k}) = strip_handles(s(j).(fn{k})); end
        end
    end
end
end
