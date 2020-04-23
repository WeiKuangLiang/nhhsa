function y = spmm_downsample(x,N,varargin)
% Shift dimension if necessary
siz = size(x);  % Save original size of x (possibly N-D)
[x,nshift] = shiftdim(x);

phase = parseUpDnSample(N,varargin{:});


    % Perform the downsample
    y = x(phase:N:end, :);

siz(nshift+1) = size(y,1);  % Update sampled dimension
y = shiftdim(y,-nshift);
if length(siz)>2
    y = reshape(y,siz);  % Restore N-D shape
end

% --------------------------------------------------------
function phase = parseUpDnSample(N,varargin)
% parseUpDnSample Parse input arguments and perform error checking.

% Initialize output args.
phase = 0;

if ( ~isnumeric(N) || (length(N) ~=1) || (fix(N) ~= N) || (N < 1) )
   error('downsample error');
end

if ~isempty(varargin)
   phase = varargin{1};
end

if ( (~isnumeric(phase)) || (fix(phase) ~= phase) || (phase > N-1) || (phase < 0))
   error('error downsample offset');
end

phase = phase + 1; % Increase phase for 1-based indexing



% [EOF]
