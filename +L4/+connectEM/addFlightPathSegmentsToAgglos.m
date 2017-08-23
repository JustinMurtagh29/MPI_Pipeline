function [ axonsNew ] = addFlightPathSegmentsToAgglos( axons, ...
    uniqueSegments, eqClassCCfull, startAgglo, idxGood )
%ADDFLIGHTPATHSEGMENTSTOAGGLOS Add single segments to agglos.
% INPUT see all inputs generated in connectEM.axonQueryAnalysis in pipeline
%           repository
%           axons: [Nx1] cell
%               Axon agglomerate equivalence classes.
%           uniqueSegments: [Nx1] cell
%               Cell array of struct containing the segments and segment
%               evidences (=occurences) of axons queries.
%           eqClassCCfull: [Nx1] cell
%               Equivalence classes of axon agglomerates.
%           startAgglo: [Nx1] cell
%               Linear index for axons where the corresponding axon query
%               started.
%           idxGood: [Nx1] logical
%               Logical indices specifying which queries were used (applies
%               to uniqueSegments, startAgglo, endAgglo).
% OUTPUT axonsNew: [Nx1] cell
%           Cell array of axon agglomerates combined by queries and with
%           single segments along flight path added. Each cell contains all
%           segment ids of one axon agglomerate.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% combine axon agglomerates from queries
axonsNew = cellfun(@(x){cell2mat(axons(x))}, eqClassCCfull);

% map from axons idx to axonsNew idx
ax2axNew = zeros(length(axons), 1, 'uint32');
ax2axNew(cell2mat(eqClassCCfull)) = ...
    repelem((1:length(eqClassCCfull)), cellfun(@length, eqClassCCfull));

% get single segments for the queries with node evidence >53
uniqueSegments = uniqueSegments(idxGood);
toAddSingleSegIds = cellfun(@(x)x.segIds(x.occurences > 53), ...
    uniqueSegments, 'uni', 0);

% get agglomerates where the single segments need to be added to
startAgglo = cell2mat(startAgglo(idxGood));
mergeIntoIdx = ax2axNew(startAgglo);
mergeIntoCC = accumarray(mergeIntoIdx, 1:length(mergeIntoIdx), [], @(x){x});
idx = mergeIntoIdx(cellfun(@(x)x(1), mergeIntoCC));

% add single segments to axonsNew
axonsNew(idx) = cellfun(@(x, y)unique([x; cell2mat(toAddSingleSegIds(y))]), ...
    axonsNew(idx), mergeIntoCC, 'uni', 0);

end

