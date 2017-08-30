function ids = getSegIds( sagglo )
%GETIDS Get the segment ids for a superagglo.
% INPUT sagglos: [Nx1] struct
%           Superagglo struct.
% OUTPUT ids: [Nx1] cell
%           Cell array of segment ids without 0.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

ids = cellfun(@(x)x(~isnan(x(:,4)), 4), {sagglo.nodes}', 'uni', 0);

end

