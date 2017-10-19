function [agglos,aggloLUT] = applyWholeCellCorrections(agglos,somaAgglos,p,folder)
% this function applied the changes (merges/splits) made in the nml files located in
% "folder" to all whole cell superagglos
%
% INPUT
% agglos       agglos in the superagglo format
% somaAgglos   agglos in the old format containing the segIds belonging to
%              each soma
% p            parameter structure from the pipeline
% folder       path to the folder containing the nml files with the changes
%
% OUTPUT
% agglos       modified agglos in the superagglo format
% aggloLUT     agglo lookup table telling which segId belongs to which
%              agglo Id
%
% by Marcel Beining <marcel.beining@brain.mpg.de>

segmentMeta = load(fullfile(p.saveFolder, 'segmentMeta.mat'),'point');

aggloSegIds = cell2mat(arrayfun(@(x) x.nodes(:,4),agglos,'uni',0));
aggloLUT(aggloSegIds)  = repelem(1:numel(agglos),arrayfun(@(x) numel(x.nodes(:,4)),agglos));

somaSegIds = cell2mat(somaAgglos);
somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
somaAgglosCoord = cellfun(@(x) segmentMeta.point(:,x)',somaAgglos,'uni',0);

% check which segId of soma is found in which agglo
[ismem,ind] = ismember(somaSegIds,aggloSegIds);
% get agglo id of each whole cell
aggloSomaId = accumarray(somaLUT(somaSegIds(ismem))',aggloLUT(aggloSegIds(ind(ismem)))',[],@mode);


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
    nodesToDelete = find(~ismember(agglos(aggloSomaId(ind)).nodes(:,1:3),cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    nodesToAdd = find(~ismember(cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)),agglos(aggloSomaId(ind)).nodes(:,1:3),'rows'));  % find node ind which has to be deleted by checking which one is missing in the loaded skeleton compared to the skeleton before modification
    
    % correct mergers
    if ~isempty(nodesToDelete)
        splitAgglo = Superagglos.removeSegIdsFromAgglos(agglos(aggloSomaId(ind)),agglos(aggloSomaId(ind)).nodes(nodesToDelete,4)); % get the splitted agglos when node is deleted
        overlapped = 0;
        segIdsToDelete = agglos(aggloSomaId(ind)).nodes(nodesToDelete,4);
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
                    agglos(end+1).nodes = [0 0 0 0];
                    aggloSomaId(n) = numel(agglos);
                end
                agglos(aggloSomaId(n)) = splitAgglo(overlapsSoma);  % replace whole cell with correctly splitted version
                aggloLUT(splitAgglo(overlapsSoma).nodes(:,4)) = repelem(aggloSomaId(n),size(splitAgglo(overlapsSoma).nodes,1));
                usedWholeCellSlots = cat(1,usedWholeCellSlots,aggloSomaId(n));
                overlapped = 1;
                overlapsSomaAll = overlapsSomaAll | overlapsSoma; % remember all split agglos that overlapped
            end
        end
        if overlapped
            % update  LUT
            aggloLUT(cell2mat(arrayfun(@(x) x.nodes(:,4),splitAgglo(~overlapsSomaAll),'uni',0))) = repelem((1:sum(~overlapsSomaAll))+numel(agglos),arrayfun(@(x) size(x.nodes,1),splitAgglo(~overlapsSomaAll)));
            agglos(end+1:end+sum(~overlapsSomaAll)) = splitAgglo(~overlapsSomaAll);  % add the splitted stuff to end of agglo class
                
            aggloLUT(segIdsToDelete) = 0;  % deleted segIds will not be there anymore in the agglo class
        else
            warning('Did not find any soma overlapping the split agglos (that should belong to soma %d) by more than 20%!',ind)
        end
    end
    
    % correct splits
    if ~isempty(nodesToAdd)
       
       skelCoords = cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0));  % putting all skel nodes together
       skelNumNodes = cellfun(@(x) size(x,1),skel.nodes);
       skelEdges = cell2mat(arrayfun(@(x) skel.edges{x}+sum(skelNumNodes(1:x-1)),(1:numel(skel.nodes))','uni',0));% putting all edges together (inceasing index for each skel
       endingSkelEdges = skelEdges(any(ismember(skelEdges,nodesToAdd),2),:);   %edges of skeletons including the last node of the original skeleton
       endingClusters = Graph.findConnectedComponents(endingSkelEdges,1,1);   % cluster each extended ending
       [~,endingSkelEdgesClusters] = cellfun(@(x) ismember(endingSkelEdges(all(ismember(endingSkelEdges,x),2),:),x),endingClusters,'uni',0);  % get the skel edges for each ending
       
%        skelSegIds = Seg.Global.getSegIds(p,skelCoords(nodesToAdd,:));  % extract the seg Ids of these nodes that were added
       % put this in later
       for n = 1:numel(endingClusters)
           skelSegIds = Seg.Global.getSegIds(p,skelCoords(endingClusters{n},:));  % extract the seg Ids of these nodes that were added
           
           indToAdd = setdiff(aggloLUT(setdiff(skelSegIds,0)),[0,aggloSomaId(ind)]); % get the index of the superagglo(s) to add
           if isempty(indToAdd)
               warning('Skel %s contained an ending which could not be processed, because the tracing did not reach a segId not belonging to the whole cell agglo')
               continue
           end
           nodesToDelete = sort(find(~ismember(skelSegIds,cell2mat(arrayfun(@(x) x.nodes(:,4),agglos([aggloSomaId(ind),indToAdd]),'uni',0)))),'descend');
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
           agglos(aggloSomaId(ind)) = Superagglos.applyEquivalences({1:numel(indToAdd)+1},agglos([aggloSomaId(ind),indToAdd]),segIdEdges);
           
%            skel2 = Superagglos.toSkel(agglos(aggloSomaId(ind)));
%            subplot(1,2,2);hold all;skel2.plot;axis equal
           aggloLUT(cell2mat(arrayfun(@(x) x.nodes(:,4),agglos(indToAdd),'uni',0))) = aggloSomaId(ind); % update LUT
           % remove agglo which has been added and update LUT
           aggloLUT = connectEM.changem(aggloLUT,(0:numel(agglos))-cumsum(accumarray(indToAdd',1,[numel(agglos)+1,1]))',0:numel(agglos));
           agglos(indToAdd) = [];
           % update aggloSomaId
           aggloSomaId =  aggloSomaId - sum(bsxfun(@ge,aggloSomaId,indToAdd),2);
       end
    end
end