function ov = getNewAxonGTAggloOverlap( gtSegIds, sagglos, ovT )
%GETNEWAXONGTAGGLOOVERLAP Get the axon agglos that overlap with the
% new_axon_gt axons.
%
% INPUT gtSegIds: [Nx1] cell
%           Segment ids of the ground truth axons. Can contain zeros and
%           duplicates.
%           (see also connectEM.eval.getNewAxonGT)
%       sagglos: [Nx1] cell or struct
%           Axon agglos or superagglos.
%       ovT: (Optional) double
%           Minimal overlap in number of segments.
%           (Default: 1)
%
% OUTPUT ov: [Nx1] cell
%           Cell with linear indices of sagglos that overlap with the
%           corresponding gtSegIds. Each cell contains a [Nx2] int array
%           where the first row is the linear index of an agglo. The second
%           row specified the number of overlapping segments.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('ovT', 'var') || isempty(ovT)
    ovT = 1;
end

% make gtSegIds unique and delete 0, -1
gtSegIds = cellfun(@(x)setdiff(x, [0, -1]), gtSegIds, 'uni', 0);

% segments ids of agglos
if isstruct(sagglos)
    sagglos = Superagglos.getSegIds(sagglos);
    sagglos = cellfun(@(x)x(~isnan(x)), sagglos, 'uni', 0);
end

% lookup table for gt seg ids in agglos
m = max([cellfun(@max, gtSegIds); ...
         cellfun(@max, sagglos(~cellfun(@isempty, sagglos)))]);
lut = Agglo.buildLUT(m, sagglos);

% get overlaps
ov = cellfun(@(x)tabulate(lut(x)), gtSegIds, 'uni', 0);
ov = cellfun(@(x)x(x(:,1) > 0, 1:2), ov, 'uni', 0); % remove id 0 as overlap
ov = cellfun(@(x)x(x(:,2) >= ovT, :), ov, 'uni', 0); % minimal overlap

end

