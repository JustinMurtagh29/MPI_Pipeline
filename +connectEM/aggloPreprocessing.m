%% load graph etc
% Comment MB: executed with clean state, e.g. ~exist will all trigger

load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');

if ~exist('graph','var') || ~all(isfield(graph,{'edges','prob','borderIdx'}))
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges','prob','borderIdx');
end
if  ~all(isfield(graph,{'neighbours','neighBorderIdx'}))
    [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
end
disp('graph loaded');

if ~exist('corrEdges','var')
    corrEdges = graph.edges(isnan(graph.borderIdx),:);
end
disp('correspondences extracted');

if ~exist('segmentMeta','var')
    segmentMeta = load(fullfile(p.saveFolder, 'segmentMeta.mat'));
end
disp('segment meta loaded');

if ~exist('borderMeta','var')
    borderMeta = load(fullfile(p.saveFolder, 'globalBorder.mat'), 'borderSize', 'borderCoM');
end
disp('border sizes loaded');

if ~exist('heuristics','var')
    heuristics = load(fullfile(p.saveFolder, 'heuristicResult.mat'));
end
disp('heuristics loaded');

%% load dendrite equivalence classes after grid search & create dendrite superagglo
if ~exist(fullfile(outputFolder,'dendrites_01.mat'),'file')
    % load all dendrite agglomerate results from after grid search
    thisGrid = load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat','axons','dendrites','heuristics','dendritesFinal','dendriteEdges');
    
    disp('Apply garbage collection')
    % Also add edges added during garbage collection (only the dendritic
    % ones), this creates the "hot edges" for the dendrite class
    [thisGrid.axonsNew, thisGrid.dendritesNew, thisGrid.garbageEdges] = connectEM.garbageCollection(graph, segmentMeta, thisGrid.axons, thisGrid.dendrites, thisGrid.heuristics.mapping);
    assert(all(cellfun(@(x,y)all(x==y), thisGrid.dendritesNew, thisGrid.dendritesFinal)));
    edgesGTall = [thisGrid.dendriteEdges; thisGrid.garbageEdges(all(ismember(thisGrid.garbageEdges,cell2mat(thisGrid.dendritesNew)),2),:)];
    save(fullfile(outputFolder,'dendritesEdgesGTall.mat'), 'edgesGTall');
    
    % get hot edges to each agglo and create first non-hybrid superagglo
    dendrites = connectEM.transformAggloOldNewRepr(thisGrid.dendritesNew,edgesGTall,segmentMeta,1);
    clear thisGrid edgesGTall;
    save(fullfile(outputFolder,'dendrites_01.mat'),'dendrites');
else    
    load(fullfile(outputFolder,'dendrites_01.mat'),'dendrites');
end

%% load/ create axon superagglos
if ~exist(fullfile(outputFolder,'axons_01.mat'),'file')
    % all axons agglomerate results from after directionality based grid
    % search (> 100 voxel)
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');
    load('/tmpscratch/mberning/edgesGTall.mat', 'edgesGTall');
    save(fullfile(outputFolder,'axonsEdgesGTall.mat'), 'edgesGTall');
    axons = connectEM.transformAggloOldNewRepr(axonsNew, edgesGTall, segmentMeta,1);
    clear axonsNew edgesGTall;
    save(fullfile(outputFolder,'axons_01.mat'),'axons');
else
    load(fullfile(outputFolder,'axons_01.mat'),'axons');
end

%% execute corresponding edges
axons = connectEM.executeEdges(axons,corrEdges,segmentMeta);
dendrites = connectEM.executeEdges(dendrites,corrEdges,segmentMeta);

save(fullfile(outputFolder,'axons_02.mat'),'axons')
save(fullfile(outputFolder,'dendrites_02.mat'),'dendrites')

%% add myelinated processes to axon class

disp('Add myelinated processes of dendrite class to axon class and execute correspondences again')
[ myelin_Dend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
myelin_Dend = myelin_Dend > 0.08;  % use the empiric myelin threshold for dendrite agglos

% add the myelinated "dendrites" to the axon class and execute
% corresponding edges between them
axons = cat(1,axons,dendrites(myelin_Dend));
dendrites = dendrites(~myelin_Dend);
fprintf('Added %d agglos of the dendritic class (now %d remaining) to the axon class (now %d agglos)',sum(myelin_Dend),numel(dendrites),numel(axons));

% execute corresponding edges again on the merged axon class. This at the
% same time merges all superagglos that have overlapping segments!
[ axons ] = connectEM.executeEdges(axons,corrEdges,segmentMeta);

%% get myelin surface scores, size scores and save final axon/dendrite class state

indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, connectEM.transformAggloNewOldRepr(dendrites));
indBigAxons = Agglo.isMaxBorderToBorderDistAbove(p, 5000, connectEM.transformAggloNewOldRepr(axons));
[ myelin_Axon ] = connectEM.calculateSurfaceMyelinScore( axons, graph, borderMeta, heuristics );  % calculate myelin score for the axon class
[ myelin_Dend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class

save(fullfile(outputFolder,'axons_03.mat'),'axons','myelin_Axon','indBigAxons')
save(fullfile(outputFolder,'dendrites_03.mat'),'dendrites','myelin_Dend','indBigDends')
