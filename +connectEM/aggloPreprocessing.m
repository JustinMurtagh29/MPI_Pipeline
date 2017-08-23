outputFolder = '/tmpscratch/mbeining/new/'  

%% load graph etc

load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat')
disp('Parameters loaded')

if ~exist('corrEdges','var')
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/correspondencesNew.mat','corrEdges');
end
disp('correspondences loaded')

if ~exist('graph','var') || ~all(isfield(graph,{'edges','prob','borderIdx'}))
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat','edges','prob','borderIdx');
end
if  ~all(isfield(graph,{'neighbours','neighBorderIdx'}))
    [graph.neighbours, graph.neighBorderIdx] = Graph.edges2Neighbors(graph.edges);
    if any(cellfun(@(x) size(x,2)>size(x,1),graph.neighbours(1:10)))  % check first 10 agglos should suffice
        graph.neighbours = cellfun(@transpose,graph.neighbours,'uni',0);
    end
end
disp('graph loaded')

if ~exist('segmentMeta','var')
    segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
end
disp('segment meta loaded')

if ~exist('borderMeta','var')
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat','borderSize', 'borderCoM');
end
disp('border sizes loaded')

if ~exist('heuristics','var')
    heuristics = load('/gaba/u/mberning/results/pipeline/20170217_ROI/heuristicResult.mat');
end
disp('heuristics loaded')

%% load / create dendrite superagglo
if ~exist(fullfile(outputFolder,'dendritesSuperPre.mat'),'file')
    thisGrid = load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat','axons','dendrites','heuristics','dendritesFinal','dendriteEdges');
    
    disp('Apply garbage collection')
    % Also add edges added during garbage collection (only the dendritic ones
    [thisGrid.axonsNew, thisGrid.dendritesNew, thisGrid.garbageEdges] = connectEM.garbageCollection(graph, segmentMeta, thisGrid.axons, thisGrid.dendrites, thisGrid.heuristics.mapping);
    assert(all(cellfun(@(x,y)all(x==y), thisGrid.dendritesNew, thisGrid.dendritesFinal)));
    edgesGTall = [thisGrid.dendriteEdges; thisGrid.garbageEdges(all(ismember(thisGrid.garbageEdges,cell2mat(thisGrid.dendritesNew)),2),:)];
    save(fullfile(outputFolder,'edgesGTall_dendrites.mat'), 'edgesGTall');
    
    % get hot edges to each agglo and create first non-hybrid superagglo
    dendrites = connectEM.transformAggloOldNewRepr(thisGrid.dendritesNew,edgesGTall,segmentMeta);
    clear thisGrid
    save(fullfile(outputFolder,'dendritesSuperPre.mat'),'dendrites')
end
%% load/ create axon superagglos
if ~exist(fullfile(outputFolder,'axonsSuperPre.mat'),'file')
    % directionality already only big axons!
    thisGrid = load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');
    thisGrid = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
    
    axons = cat(1,'axonsNew','axonsSmall');  %!!!!!!!!!!!!!!!
    axons = connectEM.transformAggloOldNewRepr(axons,graph.edges,segmentMeta);
    
    save(fullfile(outputFolder,'axonsSuperPre.mat'),'axons')
end
%% execute corresponding edges

axons = connectEM.executeEdges(axons,corrEdges);

dendrites = connectEM.executeEdges(dendrites,corrEdges);

%% add myelinated processes to axon class

disp('Add myelinated processes of dendrite class to axon class and execute correspondences again')
[ myelin_Dend ] = calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
myelin_Dend = myelin_Dend > 0.08;  % use the empiric myelin threshold for dendrite agglos

% add the myelinated "dendrites" to the axon class and execute
% corresponding edges between them
axons = cat(1,axons,dendrites(myelin_Dend));
dendrites = dendrites(~myelin_Dend);
fprintf('Added %d agglos of the dendritic class (now %d remaining) to the axon class (now %d agglos)',sum(myelin_Dend),numel(dendrites),numel(axons));
[ axons ] = connectEM.executeEdges(axons,corrEdges,segmentMeta);
[axons,n] = connectEM.removeDuplicateSegIdsInAgglo(axons);

%%

disp('Solve overlaps between added myelinated axons and axon class')
% in each agglo there is no overlap but between agglos
[ovAgglos, segId] = aggloOverlaps(axons);
[ovAgglosGrouped]  = Graph.findConnectedComponents(cell2mat(ovAgglos),1,1);

% maybe do not merge those with only one or two segment overlap..maybe also
% make a size threshold.. if doing that, remember to delete the not used 
% segments from the myelinatd axons or both

% transform the segId cell also into a grouped cell array
[~,idx] = ismember(cell2mat(ovAgglos),cat(1,ovAgglosGrouped{:}));
numAgglos = cellfun(@numel,ovAgglosGrouped);
countVec = repelem((1:numel(numAgglos))',numAgglos);  % create a indices vector to reference the idx to the ovAgglosGrouped
groupIdx = countVec(idx);  % make the idx correspond to the ovAgglosGrouped idx
assert(all(groupIdx(:,1)==groupIdx(:,2)))
segIdGrouped = arrayfun(@(x) cat(1,segId{groupIdx(:,1)==x}),1:numel(ovAgglosGrouped),'uni',0)';

ovAggloIds = unique(cat(1,ovAgglosGrouped{:}));
nonOvAggloIds = setdiff(1:numel(axons),ovAggloIds)';
equivalencesClass1 = cat(1,ovAgglosGrouped,num2cell(nonOvAggloIds));
axons = connectEM.applyEquivalences(equivalencesClass1,axons);
[axons,n] = connectEM.removeDuplicateSegIdsInAgglo(axons);

%% get myelin surface scores, size scores and save results

indBigDends = isMaxBorderToBorderDistAbove(p, 5000, connectEM.transformAggloNewOldRepr(dendrites));
indBigAxons = isMaxBorderToBorderDistAbove(p, 5000, connectEM.transformAggloNewOldRepr(axons));
[ myelin_Axon ] = calculateSurfaceMyelinScore( axons, graph, borderMeta, heuristics );  % calculate myelin score for the axon class
[ myelin_Dend ] = calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class


save(fullfile(outputFolder,'AllSuperagglosPreQuery.mat'),'dendrites','axons','myelin_Axon''myelin_Dend','indBigDends','indBigAxons')
