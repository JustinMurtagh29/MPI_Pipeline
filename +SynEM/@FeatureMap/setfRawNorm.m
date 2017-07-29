function setfRawNorm(obj, f)
% SETFRAWNORM Set the function handle for raw data normalization.
% INPUT f: function handle or string
%           The function used for raw data normalization as a string that
%           can be transformed to a function handle via str2func or as a
%           function handle that will be converted to string using
%           str2func.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isa(f, 'function_handle')
    obj.fRawNorm = func2str(f);
elseif ischar(f)
    obj.fRawNorm = f;
else
    error('Input must be a function handle or string.');
end
end
