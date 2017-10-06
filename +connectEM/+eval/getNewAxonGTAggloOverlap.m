function ov = getNewAxonGTAggloOverlap( gtSegIds, sagglos )
%GETNEWAXONGTAGGLOOVERLAP Get the axon agglos that overlap with the
% new_axon_gt axons.
% INPUT gtSegIds: [Nx1] cell
%           Segment ids of the ground truth axons. Can contain zeros and
%           duplicates.
%           (see also connectEM.eval.getNewAxonGT)
%       sagglos: [Nx1] cell or struct
%           Axon agglos or superagglos.
% OUTPUT ov: [Nx1] cell
%           Cell with linear indices of sagglos that overlap with the
%           corresponding gtSegIds.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% make gtSegIds unique and delete 0, -1
gtSegIds = cellfun(@(x)setdiff(x, [0, -1]), gtSegIds, 'uni', 0);

% segments ids of agglos
if isstruct(sagglos)
    sagglos = Superagglos.getSegIds(sagglos);
end

% lookup table for gt seg ids in agglos
m = max([cellfun(@max, gtSegIds), cellfun(@max, sagglos)]);
lut = Agglo.buildLUT(m, sagglos);

% get overlaps
ov = cellfun(@(x)unique(lut(x)), gtSegIds, 'uni', 0);

end

