function [dendrites,dendritesLUT] = applyWholeCellCorrections(dendrites,somaAgglos,p,folder,onlySplit,axons)
% this function applied the changes (merges/splits) made in the nml files located in
% "folder" to all whole cell superagglos
%
% INPUT
% dendrites       agglos in the superagglo format
% somaAgglos   agglos in the old format containing the segIds belonging to
%              each soma
% p            parameter structure from the pipeline
% folder       path to the folder containing the nml files with the changes
% onlySplit    boolean: only apply splits, not merges (DEFAULT 0)
%
% OUTPUT
% dendrites       modified agglos in the superagglo format
% dendritesLUT     agglo lookup table telling which segId belongs to which
%              agglo Id
%
% by Marcel Beining <marcel.beining@brain.mpg.de>

if ~exist('onlySplit','var') || isempty(onlySplit)
    onlySplit = 0;
end
if exist('axons','var') && ~isempty(axons)
    axonSegIds = cell2mat(arrayfun(@(x) x.nodes(:,4),axons,'uni',0));
    axonsLUT(axonSegIds)  = repelem(1:numel(axons),arrayfun(@(x) numel(x.nodes(:,4)),axons));
end


segmentMeta = load(fullfile(p.saveFolder, 'segmentMeta.mat'),'point');

dendriteSegIds = cell2mat(arrayfun(@(x) x.nodes(:,4),dendrites,'uni',0));
dendritesLUT(dendriteSegIds)  = repelem(1:numel(dendrites),arrayfun(@(x) numel(x.nodes(:,4)),dendrites));

somaSegIds = cell2mat(somaAgglos);
somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
somaAgglosCoord = cellfun(@(x) segmentMeta.point(:,x)',somaAgglos,'uni',0);

% check which segId of soma is found in which agglo
[ismem,ind] = ismember(somaSegIds,dendriteSegIds);
% get agglo id of each whole cell
aggloSomaId = accumarray(somaLUT(somaSegIds(ismem))',dendritesLUT(dendriteSegIds(ind(ismem)))',[],@mode);

if ~exist(folder,'dir')
    error('Folder %s is not existent',folder);
end
files = dir(fullfile(folder,'*.nml'));
usedCells = NaN(numel(somaAgglosCoord),1);
usedWholeCellSlots = [];
for f = 1:numel(files)
    skel = skeleton(fullfile(folder,files(f).name));  % load skeleton
    skel = skel.deleteTrees(cellfun(@numel,skel.nodes)/4==1); % delete single node trees
    % get soma index for current skeleton
    overlapsSoma = cellfun(@(y) sum(ismember(y, cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)),'rows')),somaAgglosCoord);
    [~,ind] = max(overlapsSoma);
    if ~isnan(usedCells(ind))
        error('Cell %d overlaps with skeleton of file %s but has already been used by file %s',ind,files(f).name,files(usedCells(ind)).name);
    end
    usedCells(ind) = f;
    nodesToDelete = find(~ismember(dendrites(aggloSomaId(ind)).nodes(:,1:3),cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    if ~onlySplit
        nodesToAdd = find(~ismember(cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)),dendrites(aggloSomaId(ind)).nodes(:,1:3),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    else
        nodesToAdd = [];
    end
    % correct mergers
    if ~isempty(nodesToDelete)
        splitAgglo = Superagglos.removeSegIdsFromAgglos(dendrites(aggloSomaId(ind)),dendrites(aggloSomaId(ind)).nodes(nodesToDelete,4)); % get the splitted agglos when node is deleted
        overlapped = 0;
        segIdsToDelete = dendrites(aggloSomaId(ind)).nodes(nodesToDelete,4);
        overlapsSomaAll = false(numel(splitAgglo),1);
        for n = 1:numel(somaAgglos) % go through all soma Agglos because one skeleton might have had several somata in it
            overlapsSoma = arrayfun(@(x) sum(ismember(somaAgglos{n},x.nodes(:,4))),splitAgglo)/numel(somaAgglos{n}) > 0.2;  % check which of the splitted agglos overlaps the soma by more than 20%
            if any(overlapsSoma)
                %                 skel2 = Superagglos.toSkel(splitAgglo(overlapsSoma));
                %                 subplot(1,2,2);skel2.plot;axis equal
                % this block solves the problem of a dendrite agglo with 2
                % or more merged somata in it by putting it at the back of
                % the superagglos
                if ismember(aggloSomaId(n),usedWholeCellSlots)
                    dendrites(end+1).nodes = [0 0 0 0];
                    aggloSomaId(n) = numel(dendrites);
                end
                dendrites(aggloSomaId(n)) = splitAgglo(overlapsSoma);  % replace whole cell with correctly splitted version
                dendritesLUT(splitAgglo(overlapsSoma).nodes(:,4)) = repelem(aggloSomaId(n),size(splitAgglo(overlapsSoma).nodes,1));
                usedWholeCellSlots = cat(1,usedWholeCellSlots,aggloSomaId(n));
                overlapped = 1;
                overlapsSomaAll = overlapsSomaAll | overlapsSoma; % remember all split agglos that overlapped
            end
        end
        if overlapped
            % update  LUT
            dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(:,4),splitAgglo(~overlapsSomaAll),'uni',0))) = repelem((1:sum(~overlapsSomaAll))+numel(dendrites),arrayfun(@(x) size(x.nodes,1),splitAgglo(~overlapsSomaAll)));
            dendrites(end+1:end+sum(~overlapsSomaAll)) = splitAgglo(~overlapsSomaAll);  % add the splitted stuff to end of agglo class
            
            dendritesLUT(segIdsToDelete) = 0;  % deleted segIds will not be there anymore in the agglo class
        else
            warning('Did not find any soma overlapping the split agglos (that should belong to soma %d) by more than 20%!',ind)
        end
    end
    
    % correct splits
    if ~isempty(nodesToAdd)
        [comments,treeIdx,nodeIdx] = skel.getAllComments;
        skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));  % putting all skel nodes together
        skelNumNodes = cellfun(@(x) size(x,1),skel.nodes);
        cumSkelNumNodes = cumsum(skelNumNodes);
        skelEdges = cell2mat(arrayfun(@(x) skel.edges{x}+sum(skelNumNodes(1:x-1)),(1:numel(skel.nodes))','uni',0));% putting all edges together (inceasing index for each skel
        nodeIdx = nodeIdx + cumSkelNumNodes(treeIdx);
        endingSkelEdges = skelEdges(any(ismember(skelEdges,nodesToAdd),2),:);   %edges of skeletons including the last node of the original skeleton
        endingClusters = Graph.findConnectedComponents(endingSkelEdges,1,1);   % cluster each extended ending
        [~,endingSkelEdgesClusters] = cellfun(@(x) ismember(endingSkelEdges(all(ismember(endingSkelEdges,x),2),:),x),endingClusters,'uni',0);  % get the skel edges for each ending
        
        %        skelSegIds = Seg.Global.getSegIds(p,skelCoords(nodesToAdd,:));  % extract the seg Ids of these nodes that were added
        % put this in later
        for n = 1:numel(endingClusters)
            if ~any(ismember(skelSegIds(endingSkelEdgesClusters{n}(:)),dendrites(aggloSomaId(ind)).nodes(:,4)))
                warning('Skel %s contained an ending which could not be processed, because it seemed to be part of a merged agglo which had been split away now.',skel.filename)
                continue
            end
            
            skelSegIds = Seg.Global.getSegIds(p,skelCoords(endingClusters{n},:));  % extract the seg Ids of these nodes that were added
            
            hasAxonComment = cellfun(@(x) ~isempty(strfind(x,'axon')),comments(endingClusters{n}));
            if any(hasAxonComment)
                indToAdd = setdiff(axonsLUT(setdiff(skelSegIds,0)),0); % get the index of the superagglo(s) to add
               
                if isempty(indToAdd)
                    warning('Skel %s contained an axon marked ending which could not be processed, because the tracing did not reach a segId belonging to an axon superagglo',skel.filename)
                    continue
                end
                % find nodes at segIds that are not part of the whole cell or
                % the superagglos to add and delete those
                nodesToDelete = sort(find(~ismember(skelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),cat(1,dendrites(aggloSomaId(ind)),axons(indToAdd)),'uni',0)))),'descend');
                skelSegIds(nodesToDelete) = [];
                for d = 1:numel(nodesToDelete)
                    %find neighbors
                    neighborIdx = setdiff(endingSkelEdgesClusters{n}(any(endingSkelEdgesClusters{n} == nodesToDelete(d),2),:),nodesToDelete(d));
                    endingSkelEdgesClusters{n} = unique(cat(1,endingSkelEdgesClusters{n}(~any(endingSkelEdgesClusters{n}==nodesToDelete(d),2),:),combnk(neighborIdx,2)),'rows'); % add edges bridging the deleted node
                    % delete node and reduce edge indices above this node idx by 1
                    endingClusters{n}(nodesToDelete(d)) = [];
                    endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>nodesToDelete(d)) = endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>nodesToDelete(d)) - 1;
                end
                segIdEdges = skelSegIds(endingSkelEdgesClusters{n});  % get segId edge vector of skeleton
                if size(segIdEdges,2)~=2 % fix stupid 1 value pair problem
                    segIdEdges = segIdEdges';
                end
                dendrites(aggloSomaId(ind)) = Superagglos.applyEquivalences({1:numel(indToAdd)+1},cat(1,dendrites(aggloSomaId(ind)),axons(indToAdd)),segIdEdges);
                
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(:,4),axons(indToAdd),'uni',0))) = aggloSomaId(ind); % update LUT
                % remove agglo which has been added and update LUT
                axonsLUT = connectEM.changem(axonsLUT,(0:numel(axons))-cumsum(accumarray(indToAdd',1,[numel(axons)+1,1]))',0:numel(axons));
                axons(indToAdd) = [];
                
            else                
                indToAdd = setdiff(dendritesLUT(setdiff(skelSegIds,0)),[0,aggloSomaId(ind)]); % get the index of the superagglo(s) to add
                if isempty(indToAdd)
                    warning('Skel %s contained an ending which could not be processed, because the tracing did not reach a segId not already belonging to the whole cell agglo or the segIDs were not part of the dendrite class',skel.filename)
                    continue
                end
                % find nodes at segIds that are not part of the whole cell or
                % the superagglos to add and delete those
                nodesToDelete = sort(find(~ismember(skelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),dendrites([aggloSomaId(ind),indToAdd]),'uni',0)))),'descend');
                skelSegIds(nodesToDelete) = [];
                for d = 1:numel(nodesToDelete)
                    %find neighbors
                    neighborIdx = setdiff(endingSkelEdgesClusters{n}(any(endingSkelEdgesClusters{n} == nodesToDelete(d),2),:),nodesToDelete(d));
                    endingSkelEdgesClusters{n} = unique(cat(1,endingSkelEdgesClusters{n}(~any(endingSkelEdgesClusters{n}==nodesToDelete(d),2),:),combnk(neighborIdx,2)),'rows'); % add edges bridging the deleted node
                    % delete node and reduce edge indices above this node idx by 1
                    endingClusters{n}(nodesToDelete(d)) = [];
                    endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>nodesToDelete(d)) = endingSkelEdgesClusters{n}(endingSkelEdgesClusters{n}>nodesToDelete(d)) - 1;
                end
                segIdEdges = skelSegIds(endingSkelEdgesClusters{n});  % get segId edge vector of skeleton
                if size(segIdEdges,2)~=2 % fix stupid 1 value pair problem
                    segIdEdges = segIdEdges';
                end
                dendrites(aggloSomaId(ind)) = Superagglos.applyEquivalences({1:numel(indToAdd)+1},dendrites([aggloSomaId(ind),indToAdd]),segIdEdges);
                
                %            skel2 = Superagglos.toSkel(agglos(aggloSomaId(ind)));
                %            subplot(1,2,2);hold all;skel2.plot;axis equal
                dendritesLUT(cell2mat(arrayfun(@(x) x.nodes(:,4),dendrites(indToAdd),'uni',0))) = aggloSomaId(ind); % update LUT
                % remove agglo which has been added and update LUT
                dendritesLUT = connectEM.changem(dendritesLUT,(0:numel(dendrites))-cumsum(accumarray(indToAdd',1,[numel(dendrites)+1,1]))',0:numel(dendrites));
                dendrites(indToAdd) = [];
                % update aggloSomaId
                aggloSomaId =  aggloSomaId - sum(bsxfun(@ge,aggloSomaId,indToAdd),2);
            end
        end
    end
end