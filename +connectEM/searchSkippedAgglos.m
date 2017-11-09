function searchSkippedAgglos(superagglos,outputFolder,param)
% this function searches superagglos for edges that are not part of the 
% graph which hints to skipped agglos. superagglos with such edges are
% written out as skeletons for checking
%
% marcel.beining@brain.mpg.de

if ~exist('param','var') || isempty(param)
    param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    param = param.p;
    
end

graph = load([param.saveFolder 'graphNew.mat'], 'edges');
graph.edges = sort(graph.edges,2);

for f = 1:numel(superagglos)
    % get and sort segId edges from current superagglo
    segIdEdges = reshape(superagglos(f).nodes(superagglos(f).edges,4),[],2);
%     segIdEdges(any(isnan(segIdEdges),2),:) = [];
    segIdEdges = sort(segIdEdges,2);
    
    % check which segIdEdges are not part of the graph
    ind = ~ismember(segIdEdges,graph.edges,'rows');
    if any(ind)
        % if any edge is not part of the graph there was probably one or
        % multiple agglos skipped, so make a comment there and write
        % superagglo to skelton file for check
        superagglos(f).comments = repmat({''}, size(superagglos(f).nodes,1),1);
        superagglos(f).comments(unique(superagglos(f).edges(ind,:))) = {'skipped agglo nearby'};
        skel = Superagglos.toSkel(superagglos(f));
        skel.write(fullfile(outputFolder,sprintf('Agglo_%02d.nml',f)))
    end
end