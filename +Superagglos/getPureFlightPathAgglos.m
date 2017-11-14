function idx = getPureFlightPathAgglos( sagglo )
%GETPUREFLIGHTPATHAGGLOS Get agglos that consist only of flight paths.
% INPUT sagglo: [Nx1] struct
%           Agglos in the superagglo format.
% OUTPUT idx: [Nx1] int
%           Linear indices of agglos that consist only of flight path
%           nodes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

idx = find(cellfun(@(x)all(isnan(x(:,4))), {sagglo.nodes}'));

end

