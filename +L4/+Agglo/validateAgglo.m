function [ props ] = validateAgglo( agglo )
%VALIDATEAGGLO
% INPUT agglo: [Nx1] cell
%           Cell array of intger arrays containing the ids of one
%           agglomerate per cell.
%       props: struct
%           Properties of the agglomeration:
%               isSorted
%               isUnique
%               hasExclusiveClasses
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

l = cellfun(@length, agglo);
props.isSorted = all(cellfun(@issorted, agglo));
aggloU = cellfun(@unique, agglo, 'uni', 0);
lU = cellfun(@length, aggloU);
if any(l ~= lU)
    props.isUnique = false;
else
    props.isUnique = true;
end
eClasses = Seg.Global.combineEClasses(agglo);
if length(eClasses) ~= length(agglo)
    props.hasExclusiveClasses = false;
else
    props.hasExclusiveClasses = true;
end
end

