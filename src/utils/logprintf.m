function logprintf(fid, varargin)
%LOGPRINTF  fprintf to the screen and (if fid > 0) to a log file.
fprintf(1, varargin{:});
if fid > 0, fprintf(fid, varargin{:}); end
end
