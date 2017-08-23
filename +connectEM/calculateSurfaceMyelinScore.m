function [ fracMyelinAggloBorder,aggloNeighBordersMyelin, aggloNeighBordersNoMyelin ] = calculateSurfaceMyelinScore( agglos, graph, borderMeta, heuristics )
% calculates the fraction of myelinated borders of each agglomeration
% INPUT
% agglos        either old (equivalence classes with segment Ids) or new
%               (edges and nodes) representation of superagglos
% graph         graph structure, containing at least the neighbours and
%               neighBorderIdx variable
% borderMeta    borderMeta structure containing at least the borderSize
% heuristics    heuristics of the dataset containing at least the
%               myelinScore
%
% OUTPUT
% fracMyelinAggloBorder     the fraction of myelinated borders of each agglomeration
% aggloNeighBordersMyelin   cell array comprising all myelinated border IDs of an
%                           agglo that are at the surface of the agglo
% aggloNeighBordersNoMyelin cell array comprising all non-myelinated border IDs of an
%                           agglo that are at the surface of the agglo

% transform new into old superagglo representation
agglos = connectEM.transformAggloNewOldRepr(agglos);

% get only the neighbouring segments outside of each agglo (and the indices
% in graph.neighbours
[aggloNeighbours,indaggloneighbours] = cellfun(@(x) setdiff(cat(1,graph.neighbours{x}),x),agglos,'uni',0);  

fracMyelinAggloBorder = NaN(numel(agglos),1);
aggloNeighBordersMyelin = cell(numel(agglos),1);
aggloNeighBordersNoMyelin = aggloNeighBordersMyelin;
% go through all agglos
tic
for n=1:numel(agglos)
    aggloNeighBorders = cat(1,graph.neighBorderIdx{agglos{n}}); % get all border IDs that belong to each agglomeration (outside and inbetween)
    aggloNeighBorders = aggloNeighBorders(indaggloneighbours{n});% get only those borders at the shell of each agglomeration
    ismyelin = heuristics.myelinScore(aggloNeighbours{n})>0.5;  % get index to borders which have a myelin score above 0.5
    % get indices to the borders which are no corresponding edge and either
    % myelinated or not
    aggloNeighBordersMyelin{n} = aggloNeighBorders(~isnan(aggloNeighBorders) & ismyelin);
    aggloNeighBordersNoMyelin{n} = aggloNeighBorders(~isnan(aggloNeighBorders) & ~ismyelin);
    % calculate the fraction of summed borders size that is myelinated
    fracMyelinAggloBorder(n) = sum(borderMeta.borderSize(aggloNeighBordersMyelin{n}))/sum(borderMeta.borderSize(aggloNeighBorders(~isnan(aggloNeighBorders))));
    Util.progressBar(n,numel(agglos))
end
end

