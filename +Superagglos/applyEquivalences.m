function newAgglos = applyEquivalences(equivalencesClass1,agglos,segIds)
% Apply equivalence classes on agglos (old or new representation)
% caution, due to the lookup table, every segId only has to be in a single
% agglomeration! otherwise there will be nodes without edges to it!

% INPUT
% equivalencesClass1    cell array with the equivalences of the agglos
% agglos                agglomerations (either old or new representation)
% segIds                (only for new representation) can be an edge list
%                       with segIds that belong together and for which
%                       edges are made in the superagglos
%
% OUTPUT
% newAgglos             agglos with applied equivalences
%
% author: marcel.beining@brain.mpg.de


if ~isfield(agglos,'nodes') % old representation
    newAgglos = cellfun(@(x) unique(cat(1,agglos{x})),equivalencesClass1,'uni',0);
else   % new representation, more complicated as edges have to be established
    if any(cellfun(@(x) numel(x)~=numel(unique(x)),equivalencesClass1))
        error('It seems some of the equivalence classes contain more than one reference to the same agglo. Please make sure that all Agglos contain each segment ID only once before running this function, e.g. with connectEM.removeDuplSegIdsInAgglo.')
    end
    numSegsAgglos = cellfun(@(x) size(x,1),{agglos.nodes})';
    numEdgesAgglos = cellfun(@numel,{agglos.edges})/2;
    
    newnodes = cellfun(@(x) cat(1,agglos(x).nodes),equivalencesClass1,'uni',0);
    % concatenate the edges of different agglos by adding number of
    % nodes of agglo 1:n to the edge indices of agglo n+1
    newedges = cellfun(@(x) cat(1,agglos(x).edges) + repmat(reshape(repelem(cumsum(cat(1,0,numSegsAgglos(x(1:end-1),:))),numEdgesAgglos(x)),sum(numEdgesAgglos(x)),1),1,2),equivalencesClass1,'uni',0);
    if exist('segIds','var')
        % these are normal edges with seg ID
        % create a lookup of which segmentId is found in which agglo
        lookup = zeros(max(segIds(:)),1);
        lookup(cell2mat(cellfun(@(x) x(~isnan(x(:,4)),4), newnodes, 'uni', 0))) = repelem((1:length(newnodes))', cellfun(@(x) sum(~isnan(x(:,4))), newnodes));
        aggloIdx = lookup(segIds(:,1));
        % get indices of the node to which the corr edge hints to
        [~,idx] = arrayfun(@(x) ismember(segIds(aggloIdx==x,:),newnodes{x}(:,4)),(1:numel(newnodes))','uni',0);
        % this gives the edge which has to be added to the edge list to
        % concatenate them with the corrEdge
        newedges = cellfun(@(x,y) unique(cat(1,x,sort(y(all(y,2),:),2)),'rows'),newedges,idx,'uni',0);
    end
    newAgglos = cell2struct([newedges';newnodes'],{'edges','nodes'},1);
    % remove duplicate nodes from the agglos (would be a nightmare before)
    [~,nodesToKeep] = arrayfun(@(x) unique(x.nodes,'rows'),newAgglos,'uni',0);
    newAgglos = Superagglos.removeNodesFromAgglo(newAgglos,arrayfun(@(x) setdiff(1:size(newAgglos(x).nodes,1),nodesToKeep{x}),1:numel(newAgglos),'uni',0) );
    

end

