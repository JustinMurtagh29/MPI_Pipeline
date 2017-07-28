function raw = norm2ex145( raw, mu, sig )
%NORM2EX145 Normalize raw data to mean and standard deviation of ex145.
% Normalizes the raw data input array to 122 mean and 22 std.
% INPUT raw: 3d uint8 or float
%           Raw data array.
%       mu: (Optional) float
%           Mean of raw. If you specify mu also specify sig.
%           (Default: Determined from input raw).
%       sig: (Optional) float
%           Raw standard deviation. If you specify sig also specify mu.
%           (Default: Determined from input raw).
% OUTPUT raw: 3d uint8 or float
%           Renormalized raw data.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

dtype = class(raw);
if ~exist('mu', 'var') || ~exist('sig', 'var') || isempty(mu) || isempty(sig)
    raw = zscore(single(raw));
else
    raw = (single(raw) - mu)./sig;
end
raw = raw.*22 + 122;

if strcmp(dtype, 'uint8') %#ok<STISA>
    raw = uint8(raw);
end

end
