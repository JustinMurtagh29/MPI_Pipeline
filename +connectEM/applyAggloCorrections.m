function [dendrites,dendritesLUT] = applyAggloCorrections(dendrites,p,folder,modus,axons)
% this function applied the changes (merges/splits) made in the nml files located in
% "folder" to all whole cell superagglos
%
% INPUT
% dendrites    agglos in the superagglo format
% p            parameter structure from the pipeline
% folder       path to the folder containing the nml files with the changes
% modus        0: applySplitsAndMerges; 1: only apply splits; 2: only write out agglos as skeletons that would be added by this step (to inspect for mergers) (DEFAULT 0)
% axons        agglos in the superagglo format to add axonic stuff to whole
%              cells
%
% OUTPUT
% dendrites       modified agglos in the superagglo format
% dendritesLUT     agglo lookup table telling which segId belongs to which
%              agglo Id
%
% by Marcel Beining <marcel.beining@brain.mpg.de>

if ~exist('modus','var') || isempty(modus)
    modus = 0;
end
if modus == 2
    mkdir(fullfile(folder,'checkBeforeAdd'))
    if ~isfield(dendrites,'comments')
        dendrites(1).comments = [];
    end
    if exist('axons','var') && ~isempty(axons) && ~isfield(axons,'comments')
        axons(1).comments = [];
    end
end
if exist('axons','var') && ~isempty(axons)
    axons = rmfield(axons,setdiff(fieldnames(axons),{'nodes','edges','comments'}));
    [axonsLUT] = Superagglos.buildLUT(axons);
end

[dendritesLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);
% dendriteCoord = cell2mat(arrayfun(@(x) x.nodes(:,1:3),dendrites,'uni',0));
% numDend = arrayfun(@(x) size(x.nodes,1),dendrites);
% dendLabel = repelem(1:numel(dendrites),numDend);
        
if ~exist(folder,'dir')
    error('Folder %s is not existent',folder);
end
files = dir(fullfile(folder,'*.nml'));
% usedCells = NaN(numel(somaAgglosCoord),1);
usedCells = NaN(numel(dendrites),1);

for f = 1:numel(files)
    dendritesbkp = dendrites;
    skel = skeleton(fullfile(folder,files(f).name));  % load skeleton
    skel.edges = cellfun(@(x) x(~any(x==0,2),:),skel.edges,'uni',0);

    skel = skel.deleteTrees(cellfun(@numel,skel.nodes)/4==0); % delete zero node trees, caution has been changed from single node on 08 nov 17
    skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));  % putting all skel nodes together
    % create node pairs from edges in all skels
%     skel = skel.deleteTrees(cellfun(@isempty,skel.edges)); % delete single or zero node skels
    singleNodeSkels = cellfun(@isempty,skel.edges);
    skelCoordsCatEdges = cell2mat(cellfun(@(x,y) [x(y(:,1),1:3),x(y(:,2),1:3)],skel.nodes(~singleNodeSkels),skel.edges(~singleNodeSkels),'uni',0));  % putting all skel nodes together

    warning('OFF','auxiliaryMethods:readKnossosCube')
    skelSegIds = Seg.Global.getSegIds(p,skelCoords);  % extract the seg Ids of all skel nodes
    warning('ON','auxiliaryMethods:readKnossosCube')
    % get soma index for current skeleton
%     [~,aggloOverlapsSkel] = ismember(skelCoords,dendriteCoord,'rows');
    [~,aggloOverlapsSkel] = ismember(skelSegIds,dendriteSegIds,'rows');
    ind = mode(dendritesLUT(dendriteSegIds(nonzeros(aggloOverlapsSkel))));
    
    % avoid using a wrong dendrite/axons agglo because it overlaps only a
    % little
    if sum(ismember(skelCoords,dendrites(ind).nodes(:,1:3),'rows'))/size(skelCoords,1) < 0.5
        warning('Found overlap of skeleton %s with any agglo is less than 50%%..skipping..',skel.filename);
        continue
    end
    if ~isnan(usedCells(ind))
        error('Cell %d overlaps with skeleton of file %s but has already been used by file %s',ind,files(f).name,files(usedCells(ind)).name);
    end
    usedCells(ind) = f;
    if  modus ~= 2
        % create node pairs from edges of superagglo and compare to
        % skeleton node pairs to find deleted edges
        dendCoordsCatEdges = [dendrites(ind).nodes(dendrites(ind).edges(:,1),1:3),dendrites(ind).nodes(dendrites(ind).edges(:,2),1:3)];  % putting all skel nodes together
        edgesToDelete = ~(ismember(dendCoordsCatEdges,skelCoordsCatEdges,'rows') | ismember(dendCoordsCatEdges(:,[4:6,1:3]),skelCoordsCatEdges,'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
        nodesToDelete = find(~ismember(dendrites(ind).nodes(:,1:3),skelCoords,'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    else
        nodesToDelete = [];
    end
    if  any(modus == [0 2])
        nodesToAdd = find(~ismember(skelCoords,dendrites(ind).nodes(:,1:3),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    else
        nodesToAdd = [];
    end
    % correct mergers by splitting the superagglo after edges are removed
    splitAgglo = Superagglos.removeEdgesFromAgglo(dendrites(ind),edgesToDelete); % get the splitted agglos when node is deleted
%     splitAgglo = Superagglos.removeNodesFromAgglo(dendrites(ind),nodesToDelete); % get the splitted agglos when node is deleted
    if ~isempty(nodesToDelete)
        segIdsToDelete = dendrites(ind).nodes(nodesToDelete,4);
        segIdsToDelete = segIdsToDelete(~isnan(segIdsToDelete));
    end
    dendrites(ind) = splitAgglo(1);  % replace this agglo with one of the correctly splitted version
    dendritesLUT(splitAgglo(1).nodes(~isnan(splitAgglo(1).nodes(:,4)),4)) = repelem(ind,sum(~isnan(splitAgglo(1).nodes(:,4))));
    
    % this part deletes splitted stuff that is already in the axon class
    % (to not put axonic stuff into dendrite class)
    delSplit = false(numel(splitAgglo),1);
    for s = 2:numel(splitAgglo) % go through all splitted parts of the agglo
        idxMostOVAx = mode(axonsLUT(splitAgglo(s).nodes(~isnan(splitAgglo(s).nodes(:,4)),4))); % get axon agglo idx with which the splitted part overlaps most
        if ~isnan(idxMostOVAx)
            n = sum(ismember(splitAgglo(s).nodes(:,1:3),axons(idxMostOVAx).nodes(:,1:3),'rows')); % check how much overlap incl. flight path points
            if n == size(splitAgglo(s).nodes,1) % if splitted thing fully exists in axon class, do not keep it
                delSplit(s) = true;
                if n/size(axons(idxMostOVAx).nodes,1) > 0.1
                    
                end
            end
        else  % if it is only a flight path, delete
            delSplit(s) = true;
        end
    end
    splitAgglo(delSplit) = [];
    
    
    if numel(splitAgglo) > 1
        % update  LUT
        dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),splitAgglo(2:end),'uni',0))) = repelem((1:numel(splitAgglo)-1)+numel(dendrites),arrayfun(@(x) sum(~isnan(x.nodes(:,4))),splitAgglo(2:end)));
        dendrites(end+1:end+numel(splitAgglo)-1) = splitAgglo(2:end);  % add the splitted stuff to end of agglo class
    elseif any(edgesToDelete)
        warning('Deleting the edges from the skeleton %s did not split the agglo!',skel.filename)
    end
    if ~isempty(nodesToDelete)
        dendritesLUT(segIdsToDelete) = 0;  % deleted segIds will not be there anymore in the agglo class
    end
    
    % correct splits
    if ~isempty(nodesToAdd)
        [comments,treeIdx,nodeIdx] = skel.getAllComments;
        
        skelNumNodes = cellfun(@(x) size(x,1),skel.nodes);
        cumSkelNumNodes = [0;cumsum(skelNumNodes)];
        skelEdges = cell2mat(arrayfun(@(x) skel.edges{x}+cumSkelNumNodes(x),(1:numel(skel.nodes))','uni',0));% putting all edges together (inceasing index for each skel
        nodeIdx = nodeIdx + cumSkelNumNodes(treeIdx);
        endingSkelEdges = skelEdges(any(ismember(skelEdges,nodesToAdd),2),:);   %edges of skeletons including the last node of the original skeleton
        [endingClusters,~] = Graph.findConnectedComponents(endingSkelEdges,1,1);   % cluster each extended ending
        [~,endingSkelEdgesClusters] = cellfun(@(x) ismember(endingSkelEdges(all(ismember(endingSkelEdges,x),2),:),x),endingClusters,'uni',0);  % get the skel edges for each ending
        
%         skelSegIds = Seg.Global.getSegIds(p,skelCoords(endingSkelEdges(:),:));  % extract the seg Ids of these nodes that were added
        indToAddAxonsAll = [];
        indToAddDendritesAll = [];
        for n = 1:numel(endingClusters)
%             [~,segIdInd] = ismember(endingClusters{n},endingSkelEdges(:));
            theseSkelSegIds = skelSegIds(endingClusters{n});
            if ~any(ismember(theseSkelSegIds(endingSkelEdgesClusters{n}(:)),dendrites(ind).nodes(:,4)))
                warning('Skel %s contained an ending which could not be processed, because it seemed to be part of a merged agglo which had been split away now.',skel.filename)
                continue
            end
            
            %hasAxonComment = cellfun(@(x) ~isempty(strfind(x,'axon')),comments(ismember(nodeIdx,endingClusters{n})));
            hasTrueEndingComment = cellfun(@(x) ~isempty(strfind(x,'true ending')),comments(ismember(nodeIdx,endingClusters{n})));
            if ~isempty(hasTrueEndingComment) && any(hasTrueEndingComment)
                continue
            else %~isempty(hasAxonComment) && any(hasAxonComment)
                indToAddAxons = unique(nonzeros(axonsLUT(nonzeros(theseSkelSegIds))))'; % get the index of the superagglo(s) to add
                indToAddDendrites = setdiff(dendritesLUT(nonzeros(theseSkelSegIds)),[0,ind]); % get the index of the superagglo(s) to add
                indToAddAxons = indToAddAxons(:);
                indToAddDendrites = indToAddDendrites(:);
                if all([isempty(indToAddDendrites) isempty(indToAddAxons)])
                    warning('Skel %s contained an ending which could not be processed, because the tracing did not reach a segId not already belonging to the whole cell agglo or the segIDs were not part of the dendrite/axon class',skel.filename)
                    continue
                end

                if modus == 2
                    indToAddAxonsAll = cat(1,indToAddAxonsAll,indToAddAxons);
                    for i = 1:numel(indToAddAxons)
                        [~,indComment] = ismember(theseSkelSegIds,axons(indToAddAxons(i)).nodes(:,4));
                        indComment = nonzeros(indComment);
                        if isempty(axons(indToAddAxons(i)).comments)
                            axons(indToAddAxons(i)).comments = repmat({''},size(axons(indToAddAxons(i)).nodes,1),1);
                        end
                        axons(indToAddAxons(i)).comments(indComment) = repmat({'attached segments'},numel(indComment),1);
                    end
                    indToAddDendritesAll = cat(1,indToAddDendritesAll,indToAddDendrites);
                    for i = 1:numel(indToAddDendrites)
                        [~,indComment] = ismember(theseSkelSegIds,dendrites(indToAddDendrites(i)).nodes(:,4));
                        indComment = nonzeros(indComment);
                        if isempty(dendrites(indToAddDendrites(i)).comments)
                            dendrites(indToAddDendrites(i)).comments = repmat({''},size(dendrites(indToAddDendrites(i)).nodes,1),1);
                        end
                        dendrites(indToAddDendrites(i)).comments(indComment) = repmat({'attached segments'},numel(indComment),1);
                    end
                    continue
                elseif modus == 0
                    for i = 1:numel(indToAddAxons)
                        % this part handles the problem of segment ID
                        % duplets in axon and dendrite class by either
                        % removing segIds from dendrite class that are in
                        % the axon to be added (if whole dendrite agglo is
                        % contained in the axon) or by transforming the locations
                        % in dendrite agglos with segment duplets to a flight path
                        axSegIds = axons(indToAddAxons(i)).nodes(~isnan(axons(indToAddAxons(i)).nodes(:,4)),4); % get all axon seg Ids
                        indDend = unique(dendritesLUT(axSegIds));       % get dendrites that contain these seg Ids ,too
                        if ~all(indDend==0)
                            count = histc(dendritesLUT(axSegIds),indDend);  % count how many such duplets each of the dendrite agglos have
                            if indDend(1) == 0                              % remove zeros from the LUT
                                indDend = indDend(2:end);
                                count = count(2:end);
                            end
                            canBeDeleted = arrayfun(@(x) size(x.nodes,1),dendrites(indDend))==count; % if the whole dendrite agglos is contained in the axon it can be removed from dendrite class
                            indDendRest = indDend(~canBeDeleted);  % get all dendrite agglos that have only partial overlap with the axon
                            for d = 1:numel(indDendRest) % go through these agglos, get the segId duplets and transform the axon nodes with segID duplets into a flight path
                                makeTheseNaN = ismember(axons(indToAddAxons(i)).nodes(:,4),axSegIds(indDendRest(d) == dendritesLUT(axSegIds)));
                                axonsLUT(axons(indToAddAxons(i)).nodes(makeTheseNaN,4)) = 0; % remove ref to the axon for these seg ids
                                axons(indToAddAxons(i)).nodes(makeTheseNaN,4) = NaN;
                                axons(indToAddAxons(i)).nodes(makeTheseNaN,1:3) = axons(indToAddAxons(i)).nodes(makeTheseNaN,1:3)+0.1; % add tiny value to coordinate to make it different from segment centroid
                            end
                            if any(canBeDeleted)
                                % update indices to agglomerate
                                ind = ind - sum(indDend(canBeDeleted) <= ind);
                                indToAddDendrites = setdiff(indToAddDendrites,indDend(canBeDeleted));
                                indToAddDendrites = connectEM.changem(indToAddDendrites,((1:numel(dendrites))-cumsum(accumarray(indDend(canBeDeleted),1,[numel(dendrites),1]))'),1:numel(dendrites));
                                [dendrites,dendritesLUT] = Superagglos.remove(dendrites,indDend(canBeDeleted),dendritesLUT);
                            end
                        end
                    end
                    % recheck which axons will be added, as some might get
                    % thrown out now 
                    indToAddAxons = unique(nonzeros(axonsLUT(nonzeros(theseSkelSegIds))))'; % get the index of the superagglo(s) to add
                    indToAddAxons = indToAddAxons(:);
                end

                % find nodes at segIds that are not part of the whole cell or
                % the superagglos to add and delete those
                endingNodesToDelete = sort(find(~ismember(theseSkelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),cat(1,dendrites([ind;indToAddDendrites]),axons(indToAddAxons)),'uni',0)))),'descend');
                theseSkelSegIds(endingNodesToDelete) = [];
                for d = 1:numel(endingNodesToDelete)
                    %find neighbors
                    neighborIdx = setdiff(endingSkelEdgesClusters{n}(any(endingSkelEdgesClusters{n} == endingNodesToDelete(d),2),:),endingNodesToDelete(d));
                    endingSkelEdgesClusters{n} = unique(cat(1,endingSkelEdgesClusters{n}(~any(endingSkelEdgesClusters{n}==endingNodesToDelete(d),2),:),combnk(neighborIdx,2)),'rows'); % add edges bridging the deleted node
                    % delete node and reduce edge indices above this node idx by 1
                    endingClusters{n}(endingNodesToDelete(d)) = [];
                    endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>endingNodesToDelete(d)) = endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>endingNodesToDelete(d)) - 1;
                end

                segIdEdges = theseSkelSegIds(endingSkelEdgesClusters{n});  % get segId edge vector of skeleton
                if size(segIdEdges,2)~=2 % fix stupid 1 value pair problem
                    segIdEdges = segIdEdges';
                end
                dendrites(ind) = Superagglos.applyEquivalences({1:numel(indToAddAxons)+numel(indToAddDendrites)+1},cat(1,dendrites([ind;indToAddDendrites]),axons(indToAddAxons)),segIdEdges);
                
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),axons(indToAddAxons),'uni',0))) = ind; % update LUT
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),dendrites(indToAddDendrites),'uni',0))) = ind; % update LUT
                 % remove agglo which has been added and update LUT
                axonsLUT = connectEM.changem(axonsLUT,(0:numel(axons))- [0, cumsum(accumarray(indToAddAxons,1,[numel(axons),1]))'],0:numel(axons));
                axons(indToAddAxons) = [];
                [dendrites,dendritesLUT] = Superagglos.remove(dendrites,indToAddDendrites,dendritesLUT);
                ind = ind - sum(indToAddDendrites <= ind); % update index to agglomerate
                
            end
        end
%         segIds = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
%         assert(numel(segIds)==numel(unique(segIds)))
        if modus == 2
            % write out all axons and skeletons that would be added this
            % turn
            if ~isempty(indToAddAxonsAll)
                for i = 1:numel(indToAddAxonsAll)
                    connectEM.generateSkeletonFromAggloNew(axons(indToAddAxonsAll(i)),{sprintf('AxonsToBeAdded_%d',indToAddAxonsAll(i))} , fullfile(folder,'checkBeforeAdd'), [],[],sprintf('AxonsToBeAdded_%d.nml',indToAddAxonsAll(i)));
                end
            end
            if ~isempty(indToAddDendritesAll)
                for i = 1:numel(indToAddDendritesAll)
                    connectEM.generateSkeletonFromAggloNew(dendrites(indToAddDendritesAll(i)),{sprintf('DendritesToBeAdded_%d',indToAddDendritesAll(i))} , fullfile(folder,'checkBeforeAdd'), [],[],sprintf('DendritesToBeAdded_%d.nml',indToAddDendritesAll(i)));
                end
            end
        end
    end
end