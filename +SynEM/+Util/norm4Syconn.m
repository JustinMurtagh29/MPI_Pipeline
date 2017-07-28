function raw = norm4Syconn(raw, mu, sig)
%NORM4SYCONN Raw data normalization to (0, 1) trained on the ex145 data.
% INPUT raw: 3d uint8 or float
%           3d raw data matrix.
%       mu: (Optional) float
%           Mean of raw. If you specify mu also specify sig.
%           (Default: determined from data)
%       sig: (Optional) float
%           Stdandard deviation of raw. If you specify sig also specify mu.
%           (Default: determined from data)
% OUTPUT raw: 3d float
%           Data normalized to (0, 1).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mu', 'var') || ~exist('sig', 'var') || isempty(mu) || isempty(sig)
    mu = [];
    sig = [];
end

raw = SynEM.Util.norm2ex145(single(raw), mu, sig);
raw = raw./255;
end
