function [dendrites,dendritesLUT] = applyAggloCorrections(dendrites,p,folder,modus,axons)
% this function applied the changes (merges/splits) made in the nml files located in
% "folder" to all whole cell superagglos
%
% INPUT
% dendrites    agglos in the superagglo format
% p            parameter structure from the pipeline
% folder       path to the folder containing the nml files with the changes
% mode         0: applySplitsAndMerges; 1: only apply splits; 2: only write out agglos as skeletons that would be added by this step (to inspect for mergers) (DEFAULT 0)
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
    axons = rmfield(axons,'endings');
    axonSegIds = cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),axons,'uni',0));
    axonsLUT(axonSegIds(axonSegIds))  = repelem(1:numel(axons),arrayfun(@(x) numel(x.nodes(~isnan(x.nodes(:,4)),4)),axons));
end


dendriteSegIds = cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),dendrites,'uni',0));
dendritesLUT(dendriteSegIds)  = repelem(1:numel(dendrites),arrayfun(@(x) numel(x.nodes(~isnan(x.nodes(:,4)),4)),dendrites));
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

    skel = skeleton(fullfile(folder,files(f).name));  % load skeleton
    skel = skel.deleteTrees(cellfun(@numel,skel.nodes)/4==1); % delete single node trees
    skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));  % putting all skel nodes together
    skelSegIds = Seg.Global.getSegIds(p,skelCoords);  % extract the seg Ids of all skel nodes
    
    % get soma index for current skeleton
%     [~,aggloOverlapsSkel] = ismember(skelCoords,dendriteCoord,'rows');
    [~,aggloOverlapsSkel] = ismember(skelSegIds,dendriteSegIds,'rows');
    aggloOverlapsSkel = setdiff(aggloOverlapsSkel,0);
    ind = mode(dendritesLUT(dendriteSegIds(aggloOverlapsSkel)));
    
    % avoid using a wrong dendrite/axons agglo because it overlaps only a
    % little
    if sum(dendritesLUT(dendriteSegIds(aggloOverlapsSkel))==ind)/size(skelCoords,1) < 0.5
        warning('Found overlap of skeleton %s with an agglo is less than 50%..skipping..',skel.filename);
        continue
    end
    if ~isnan(usedCells(ind))
        error('Cell %d overlaps with skeleton of file %s but has already been used by file %s',ind,files(f).name,files(usedCells(ind)).name);
    end
    usedCells(ind) = f;
    if  modus ~= 2
        % refresh if there were splits done
%         dendriteCoord = cell2mat(arrayfun(@(x) x.nodes(:,1:3),dendrites,'uni',0));
%         numDend = arrayfun(@(x) size(x.nodes,1),dendrites);
%         dendLabel = repelem(1:numel(dendrites),numDend);
        nodesToDelete = find(~ismember(dendrites(ind).nodes(:,1:3),skelCoords,'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    else
        nodesToDelete = [];
    end
    if  any(modus == [0 2])
        nodesToAdd = find(~ismember(skelCoords,dendrites(ind).nodes(:,1:3),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    else
        nodesToAdd = [];
    end
    % correct mergers
    if ~isempty(nodesToDelete)
        splitAgglo = Superagglos.removeNodesFromAgglo(dendrites(ind),nodesToDelete); % get the splitted agglos when node is deleted
        segIdsToDelete = dendrites(ind).nodes(nodesToDelete,4);
      
        dendrites(ind) = splitAgglo(1);  % replace this agglo with one of the correctly splitted version
        dendritesLUT(splitAgglo(1).nodes(~isnan(splitAgglo(1).nodes(:,4)),4)) = repelem(ind,sum(~isnan(splitAgglo(1).nodes(:,4))));

        if numel(splitAgglo) > 1
            % update  LUT
            dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),splitAgglo(2:end),'uni',0))) = repelem((1:numel(splitAgglo)-1)+numel(dendrites),arrayfun(@(x) sum(~isnan(x.nodes(:,4))),splitAgglo(2:end)));
            dendrites(end+1:end+numel(splitAgglo)-1) = splitAgglo(2:end);  % add the splitted stuff to end of agglo class
            dendritesLUT(segIdsToDelete) = 0;  % deleted segIds will not be there anymore in the agglo class
        else
            warning('Deleting the nodes from the skeleton %s did not split the agglo!',skel.filename)
        end
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
        
        indToAddAxons = [];
        indToAddDendrites = [];
        for n = 1:numel(endingClusters)
%             [~,segIdInd] = ismember(endingClusters{n},endingSkelEdges(:));
            theseSkelSegIds = skelSegIds(endingClusters{n});
            if ~any(ismember(theseSkelSegIds(endingSkelEdgesClusters{n}(:)),dendrites(ind).nodes(:,4)))
                warning('Skel %s contained an ending which could not be processed, because it seemed to be part of a merged agglo which had been split away now.',skel.filename)
                continue
            end
            
            hasAxonComment = cellfun(@(x) ~isempty(strfind(x,'axon')),comments(ismember(nodeIdx,endingClusters{n})));
            
            hasTrueEndingComment = cellfun(@(x) ~isempty(strfind(x,'true ending')),comments(ismember(nodeIdx,endingClusters{n})));
            if ~isempty(hasTrueEndingComment) && any(hasTrueEndingComment)
                continue
            elseif ~isempty(hasAxonComment) && any(hasAxonComment)
                indToAdd = setdiff(axonsLUT(setdiff(theseSkelSegIds,0)),0); % get the index of the superagglo(s) to add
               
                if isempty(indToAdd)
                    warning('Skel %s contained an axon marked ending which could not be processed, because the tracing did not reach a segId belonging to an axon superagglo',skel.filename)
                    continue
                end
                
                if modus == 2
                   indToAddAxons = cat(2,indToAddAxons,indToAdd);
                   for i = 1:numel(indToAdd)
                       [~,indComment] = ismember(theseSkelSegIds,axons(indToAdd(i)).nodes(:,4));
                       indComment = setdiff(indComment,0);
                       if isempty(axons(indToAdd(i)).comments)
                           axons(indToAdd(i)).comments = repmat({''},size(axons(indToAdd(i)).nodes,1),1);
                       end
                       axons(indToAdd(i)).comments(indComment) = repmat({'attached segments'},numel(indComment),1);
                   end
                   continue 
                end
                % find nodes at segIds that are not part of the whole cell or
                % the superagglos to add and delete those
                endingNodesToDelete = sort(find(~ismember(theseSkelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),cat(1,dendrites(ind),axons(indToAdd)),'uni',0)))),'descend');
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
                dendrites(ind) = Superagglos.applyEquivalences({1:numel(indToAdd)+1},cat(1,dendrites(ind),axons(indToAdd)),segIdEdges);
                
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),axons(indToAdd),'uni',0))) = ind; % update LUT
                % remove agglo which has been added and update LUT
                axonsLUT = connectEM.changem(axonsLUT,(0:numel(axons))-cumsum(accumarray(indToAdd',1,[numel(axons)+1,1]))',0:numel(axons));
                axons(indToAdd) = [];
                
            else                
                indToAdd = setdiff(dendritesLUT(setdiff(theseSkelSegIds,0)),[0,ind]); % get the index of the superagglo(s) to add
                if isempty(indToAdd)
                    warning('Skel %s contained an ending which could not be processed, because the tracing did not reach a segId not already belonging to the whole cell agglo or the segIDs were not part of the dendrite class',skel.filename)
                    continue
                end
                
                if modus == 2
                   indToAddDendrites = cat(2,indToAddDendrites,indToAdd);
                   for i = 1:numel(indToAdd)
                       [~,indComment] = ismember(theseSkelSegIds,dendrites(indToAdd(i)).nodes(:,4));
                       indComment = setdiff(indComment,0);
                       if isempty(dendrites(indToAdd(i)).comments)
                           dendrites(indToAdd(i)).comments = repmat({''},size(dendrites(indToAdd(i)).nodes,1),1);
                       end
                       dendrites(indToAdd(i)).comments(indComment) = repmat({'attached segments'},numel(indComment),1);
                   end
                   continue 
                end
                % find nodes at segIds that are not part of the whole cell or
                % the superagglos to add and delete those
                endingNodesToDelete = sort(find(~ismember(theseSkelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),dendrites([ind,indToAdd]),'uni',0)))),'descend');
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
                dendrites(ind) = Superagglos.applyEquivalences({1:numel(indToAdd)+1},dendrites([ind,indToAdd]),segIdEdges);
                
                %            skel2 = Superagglos.toSkel(agglos(ind));
                %            subplot(1,2,2);hold all;skel2.plot;axis equal
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(~isnan(x.nodes(:,4)),4),dendrites(indToAdd),'uni',0))) = ind; % update LUT
                % remove agglo which has been added and update LUT
                dendritesLUT = connectEM.changem(dendritesLUT,(0:numel(dendrites))-cumsum(accumarray(indToAdd',1,[numel(dendrites)+1,1]))',0:numel(dendrites));
                dendrites(indToAdd) = [];
            end
        end
        if modus == 2
            % write out all axons and skeletons that would be added this
            % turn
            connectEM.generateSkeletonFromAggloNew(cat(1,dendrites(indToAddDendrites),axons(indToAddAxons)),arrayfun(@(x) sprintf('AggloToBeAdded_%02d',x),1:numel(cat(2,indToAddDendrites,indToAddAxons)),'uni',0) , fullfile(folder,'checkBeforeAdd'), [],[],sprintf('checkBeforeAdd_%02d.nml',f));
        end
    end
end