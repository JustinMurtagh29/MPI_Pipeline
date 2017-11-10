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
graph.edges = unique(graph.edges,'rows');

% create adjacency matrix (= LUT)
maxSeg = double(max(graph.edges(:)));
am = sparse( ...
    double(graph.edges(:, 1)'), double(graph.edges(:, 2)'), ...
    true(size(graph.edges,1),1), maxSeg, maxSeg);

% make symmetric
am = am + am';

    
for f = 1:numel(superagglos)
    % get and sort segId edges from current superagglo
    segIdEdges = reshape(superagglos(f).nodes(superagglos(f).edges,4),[],2);
    segIdEdges(any(isnan(segIdEdges),2),:) = [];
%     segIdEdges = sort(segIdEdges,2);
    
    % check which segIdEdges are not part of the graph
    ind = find(am(sub2ind(size(am),segIdEdges(:,1),segIdEdges(:,2)))==0);

%     ind = find(~ismember(segIdEdges,graph.edges,'rows'));
    if ~isempty(ind)
        % delete found edges and check if there is no other path between
        % the two nodes. If this is the case it is a real occasion that
        % should be checked
        skel = Superagglos.toSkel(superagglos(f));
        skel.edges{1}(ind,:) = [];
        indReal = [];
        for e = 1:numel(ind)
            if isempty(skel.getShortestPath(superagglos(f).edges(ind(e),1), superagglos(f).edges(ind(e),2)))
                indReal = cat(1,indReal,ind(e));
            end
        end
        % if any edge is not part of the graph there was probably one or
        % multiple agglos skipped, so make a comment there and write
        % superagglo to skelton file for check
        superagglos(f).comments = repmat({''}, size(superagglos(f).nodes,1),1);
        % only take one of the two edges and add a comment there
        superagglos(f).comments(unique(superagglos(f).edges(indReal,1))) = {'skipped agglo nearby'};
        skel = Superagglos.toSkel(superagglos(f));
        skel.write(fullfile(outputFolder,sprintf('Agglo_%02d.nml',f)))
    end
end