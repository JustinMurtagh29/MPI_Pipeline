function [ ovAgglos, segId ] = aggloOverlaps( agglo )
%AGGLOOVERLAPS Find agglomerates that overlap.
% INPUT agglo: [Nx1] cell
%           Cell array of integer arrays containing the ids of one
%           agglomerate per cell.
% OUTPUT ovAgglos: [Nx1] cell
%           Cell array containing the linear indices of overlapping
%           agglomerations in each cell.
%        segIds: [Nx1] cell
%           Segment ids that overlap in the corresponding ovAgglos.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

agglo = agglo(:);
l = cellfun(@length, agglo);
aggloId = repelem(1:length(agglo), l)';
id = cell2mat(agglo);
isInAggloId = accumarray(id, aggloId, [], @(x){x});
isInAggloId = isInAggloId(id);
l = cellfun(@length, isInAggloId);
ovAgglos = isInAggloId(l > 1);
ovAgglos = cellfun(@sort, ovAgglos, 'uni', 0);
segId = id(l > 1);


% get unique overlaps
l = cellfun(@length, ovAgglos);
outCell = cell(length(ovAgglos), 1);
outSegIds = cell(length(ovAgglos), 1);
count = 1;
for i = 2:max(l)
    curOVs = cell2mat(ovAgglos(l == i)')';
    [curOvs, ~, ic] = unique(curOVs, 'rows');
    curSegIds = accumarray(ic, segId(l == i), [], @(x){x});
    curSegIds = cellfun(@unique, curSegIds, 'uni', 0);
    outCell(count:count + size(curOvs, 1) - 1) = num2cell(curOvs, 2);
    outSegIds(count:count + size(curOvs, 1) - 1) = curSegIds;
    count = count + length(curOvs);
end
ovAgglos = outCell(1:count-1);
segId = outSegIds(1:count-1);

end

