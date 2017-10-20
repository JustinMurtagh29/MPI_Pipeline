%% load graph etc
% Comment MB: executed with clean state, e.g. ~exist will all trigger
overwrite = 0;
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');

if ~exist('graph','var') || ~all(isfield(graph,{'edges','prob','borderIdx'}))
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges','prob','borderIdx');
end
if  ~all(isfield(graph,{'neighbours','neighBorderIdx'}))
    [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
    clear neighboursIdx
end
disp('graph loaded');

if ~exist('corrEdges','var')
    corrEdges = graph.edges(isnan(graph.borderIdx),:);
end
disp('correspondences extracted');

if ~exist('segmentMeta','var')
    segmentMeta = load(fullfile(p.saveFolder, 'segmentMeta.mat'));
end
if ~isfield(segmentMeta,'dendriteProb')
    segmentMeta = connectEM.addSegmentClassInformation(p,segmentMeta);
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

% remove myelinated segments from corrEdges
corrEdges = corrEdges(all(heuristics.myelinScore(corrEdges) <= 0.5,2),:);

%% load dendrite equivalence classes after grid search & create dendrite superagglo
if ~exist(fullfile(outputFolder,'dendrites_01.mat'),'file') || overwrite
    % load all dendrite agglomerate results from after grid search
    thisGrid = load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat','axons','dendrites','heuristics','dendritesFinal','dendriteEdges');
    
    disp('Apply garbage collection')
    % Also add edges added during garbage collection (only the dendritic
    % ones), this creates the "hot edges" for the dendrite class
    [thisGrid.axonsNew, thisGrid.dendritesNew, thisGrid.garbageEdges] = connectEM.garbageCollection(graph, segmentMeta, thisGrid.axons, thisGrid.dendrites, thisGrid.heuristics.mapping);
    assert(all(cellfun(@(x,y)all(x==y), thisGrid.dendritesNew, thisGrid.dendritesFinal)));
    edgesGTall = [thisGrid.dendriteEdges; thisGrid.garbageEdges(all(ismember(thisGrid.garbageEdges,cell2mat(thisGrid.dendritesNew)),2),:)];
    % remove all duplicates from hot edge list
    edgesGTall = unique(sort(edgesGTall,2),'rows');
    save(fullfile(outputFolder,'dendritesEdgesGTall.mat'), 'edgesGTall');
    
    % add the single segments
    singleSegDendrites = (double(setdiff(segmentMeta.segIds(segmentMeta.dendriteProb > 0.3),cell2mat(thisGrid.dendritesNew))));
    
    % concatenate agglos and single segments but first throw out all
    % myelinated segments that might have been added by garbage collection
    % etc
    dendritesNoMyelin = cat(1, cellfun(@(x) x(heuristics.myelinScore(x) <= 0.5),thisGrid.dendritesNew,'uni',0),num2cell(singleSegDendrites(heuristics.myelinScore(singleSegDendrites) <= 0.5)));
    dendritesNoMyelin = dendritesNoMyelin(~cellfun(@isempty,dendritesNoMyelin));
    % recheck agglos as there might be unconnected parts in some agglos
    dendVec = cell2mat(dendritesNoMyelin);
    edgesFiltered = edgesGTall(all(ismember(edgesGTall,dendVec),2),:);
    dendrites = Graph.findConnectedComponents(cat(1,edgesFiltered,repmat(dendVec,1,2)),0,1);
    
    % get hot edges to each agglo and create first non-hybrid superagglo
    dendrites = Superagglos.transformAggloOldNewRepr(dendrites,edgesGTall,segmentMeta,1);
    indSingleSeg = true(numel(dendrites),1);
    indSingleSeg(1:end-numel(singleSegDendrites)) = false;
    clear thisGrid edgesGTall dendritesNoMyelin edgesFiltered dendVec;
    save(fullfile(outputFolder,'dendrites_01.mat'),'dendrites','indSingleSeg');
else
    load(fullfile(outputFolder,'dendrites_01.mat'),'dendrites');
end

%% load/ create axon superagglos
if ~exist(fullfile(outputFolder,'axons_01.mat'),'file') || overwrite
    % all axons agglomerate results from after directionality based grid
    % search (> 100 voxel)
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');
    load('/tmpscratch/mberning/edgesGTall.mat', 'edgesGTall');
    % remove all duplicates from hot edge list
    edgesGTall = unique(sort(edgesGTall,2),'rows');
    save(fullfile(outputFolder,'axonsEdgesGTall.mat'), 'edgesGTall');
    % throw out all myelinated segments that might have been added in
    %garbage collection etc
    axonsNoMyelin = cellfun(@(x) x(heuristics.myelinScore(x) <= 0.5),axonsNew,'uni',0);
    axonsNoMyelin = axonsNoMyelin(~cellfun(@isempty,axonsNoMyelin));
    % recheck agglos as there might be unconnected parts in some agglos
    axVec = cell2mat(axonsNoMyelin);
    edgesFiltered = edgesGTall(all(ismember(edgesGTall,axVec),2),:);
    axons = Graph.findConnectedComponents(cat(1,edgesFiltered,repmat(axVec,1,2)),0,1);
    axons = Superagglos.transformAggloOldNewRepr(axons, edgesGTall, segmentMeta,1);
    clear axonsNew edgesGTall axonsNoMyelin axVec edgesFiltered;
    save(fullfile(outputFolder,'axons_01.mat'),'axons');
else
    load(fullfile(outputFolder,'axons_01.mat'),'axons');
end
disp('State 01 superagglos loaded/generated')

%% execute corresponding edges
if ~exist(fullfile(outputFolder,'axons_02.mat'),'file') || overwrite
    axons = connectEM.executeEdges(axons,corrEdges,segmentMeta);
    save(fullfile(outputFolder,'axons_02.mat'),'axons');
else
    load(fullfile(outputFolder,'axons_02.mat'),'axons');
end
if ~exist(fullfile(outputFolder,'dendrites_02.mat'),'file') || overwrite
    dendrites = connectEM.executeEdges(dendrites,corrEdges,segmentMeta);
    indSingleSeg = arrayfun(@(x) size(x.nodes,1),dendrites)==1;
    save(fullfile(outputFolder,'dendrites_02.mat'),'dendrites','indSingleSeg');
else
    load(fullfile(outputFolder,'dendrites_02.mat'),'dendrites','indSingleSeg');
end
disp('State 02 superagglos loaded/generated')

% %% remove somata from all classes and put into separate class
% if ~exist(fullfile(outputFolder,'somata_02b.mat'),'file')
%     load('/gaba/u/rhesse/forBenedikt/somasNoClosedHoles.mat')
%     ids = cell2mat(somas(:,3));
%     dendrites = Superagglos.removeSegIdsFromAgglos(dendrites,ids);
%     axons = Superagglos.removeSegIdsFromAgglos(axons,ids);
%     somata = Superagglos.transformAggloOldNewRepr(somas(:,3), graph.edges, segmentMeta,1);
%     
%     save(fullfile(outputFolder,'axons_02b.mat'),'axons');
%     save(fullfile(outputFolder,'dendrites_02b.mat'),'dendrites');
%     save(fullfile(outputFolder,'somata_02b.mat'),'somata');
% else
%     load(fullfile(outputFolder,'axons_02b.mat'),'axons');
%     load(fullfile(outputFolder,'dendrites_02b.mat'),'dendrites');
%     load(fullfile(outputFolder,'somata_02b.mat'),'somata');
% end

%% add myelinated processes to axon class
if ~exist(fullfile(outputFolder,'axons_03.mat'),'file') || ~exist(fullfile(outputFolder,'dendrites_03_v2.mat'),'file') || overwrite
    disp('Add myelinated processes of dendrite class to axon class and execute correspondences again')
    
    % load somata including merged somata
    somaAgglos = load(fullfile(outputFolder,'somas_with_merged_somas.mat'));
    somaAggloIds = cell2mat(somaAgglos.somas(:,3));
    
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class    
%     [dendriteProbDend,axonProbDend,gliaProbDend,voxSize,hasSoma] = arrayfun(@(x) deal(median(segmentMeta.dendriteProb(x.nodes(:,4))),median(segmentMeta.axonProb(x.nodes(:,4))),median(segmentMeta.gliaProb(x.nodes(:,4))),sum(segmentMeta.voxelCount(x.nodes(:,4))),any(ismember(x.nodes(:,4),somaAggloIds))),dendrites);
    [numSeg,voxSize] = arrayfun(@(x) deal(size(x.nodes,1),sum(segmentMeta.voxelCount(x.nodes(:,4)))),dendrites);
    dendLUT = repelem(1:numel(dendrites),numSeg);
    hasSoma = false(numel(dendrites),1);
    hasSoma(dendLUT(ismember(cell2mat(Superagglos.transformAggloNewOldRepr(dendrites)),somaAggloIds))) = true;
    % use the empiric myelin threshold for dendrite agglos and other
    % thresholds to get the myelinated axons
    moveToAxon = ~hasSoma & ((numSeg > 25 & myelinDend > 0.08) | (myelinDend > 0.25)) & voxSize > 200000;% & axonProbDend >= dendriteProbDend 
    indMoveToAxon = find(moveToAxon);
    connectEM.generateSkeletonFromAggloNew(dendrites(moveToAxon), arrayfun(@(x) sprintf('dendSkel_%d',indMoveToAxon(x)),1:numel(indMoveToAxon),'uni',0), '/tmpscratch/mbeining/myelinMove/', segmentMeta.maxSegId)


    % corresponding edges between them
    axons = cat(1,axons,dendrites(moveToAxon));
    dendrites = dendrites(~moveToAxon);
    fprintf('Added %d agglos of the dendritic class (now %d remaining) to the axon class (now %d agglos)\n',sum(moveToAxon),numel(dendrites),numel(axons));
    
    % execute corresponding edges again on the merged axon class. This at the
    % same time merges all superagglos that have overlapping segments!
    [ axons ] = connectEM.executeEdges(axons,corrEdges,segmentMeta);
    
    disp('Overlaps solved, now finding 5 micron agglos and surface myelin scores');
    % get myelin surface scores, size scores and save final axon/dendrite class state
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    indBigAxons = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(axons));
    [ myelinAxon ] = connectEM.calculateSurfaceMyelinScore( axons, graph, borderMeta, heuristics );  % calculate myelin score for the axon class
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    
    save(fullfile(outputFolder,'axons_03.mat'),'axons','myelinAxon','indBigAxons');
    save(fullfile(outputFolder,'dendrites_03_v2.mat'),'dendrites','myelinDend','indBigDends');
else
    load(fullfile(outputFolder,'axons_03.mat'),'axons','myelinAxon','indBigAxons');
    load(fullfile(outputFolder,'dendrites_03_v2.mat'),'dendrites','myelinDend','indBigDends');
end
disp('State 03 superagglos loaded/generated')

%% load soma whole cell agglos
somaAgglos = connectEM.getCenterSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'));
somaAgglos = somaAgglos(~cellfun(@isempty,somaAgglos)); % remove not found somata
somaSegIds = cell2mat(somaAgglos);
% remove duplicate segIds
[~,ic] = unique(somaSegIds);
duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
% somaDuplicateIds = cellfun(@(x) any(intersect(x,duplicates)),somaAgglos);
somaAgglos = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
somaSegIds = cell2mat(somaAgglos);
somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
disp('Soma whole cell agglos loaded')


%% apply manual fixes of whole cells in dendrite class
if ~exist(fullfile(outputFolder,'dendrites_03_v2_splitmerged.mat'),'file') || overwrite
    correctionFolder = 'WholeCellCorrections_03_v2';
    fprintf('Folder with correction nmls for state dendrites_03_v2 is %s\n',fullfile(outputFolder,correctionFolder));
%     [dendrites,dendriteLUT] = connectEM.applyWholeCellCorrections(dendrites,somaAgglos,p,fullfile(outputFolder,correctionFolder),1);
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),1);
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_03_v2_splitmerged.mat'),'dendrites','WholeCellId','myelinDend','indBigDends');
else
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_03_v2_splitmerged.mat'),'dendrites','WholeCellId')%,'myelinDend','indBigDends');
end
disp('Dendrites state 03 splitmerged superagglos loaded/generated')

%% next round of fixes
% apply manual fixes of whole cells in dendrite class
if ~exist(fullfile(outputFolder,'dendrites_04.mat'),'file') || overwrite
    correctionFolder = 'WholeCellCorrections_03_v2_splitmerged';
    fprintf('Folder with correction nmls for state dendrites_03_v2_splitmerged is %s\n',fullfile(outputFolder,correctionFolder));
%     [dendrites,dendriteLUT] = connectEM.applyWholeCellCorrections(dendrites,somaAgglos,p,fullfile(outputFolder,correctionFolder));
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder));

    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_04.mat'),'dendrites','WholeCellId','myelinDend','indBigDends');
else
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_04.mat'))
end
disp('Dendrites state 04 superagglos loaded/generated')
%% next round of fixes
% apply manual fixes of whole cells in dendrite class
if ~exist(fullfile(outputFolder,'dendrites_05.mat'),'file') || overwrite
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    correctionFolder = 'WholeCellCorrections_04';
    fprintf('Folder with correction nmls for state dendrites_04 is %s\n',fullfile(outputFolder,correctionFolder));
%     [dendrites,dendriteLUT] = connectEM.applyWholeCellCorrections(dendrites,somaAgglos,p,fullfile(outputFolder,correctionFolder),axons);
    dendrites = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now
    dendrites = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder,'checkBeforeAdd'),1);
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkBeforeAdd'),1);
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    wholeCells = dendrites(WholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),WholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_01.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_05.mat'),'dendrites','WholeCellId','myelinDend','indBigDends');
    
else
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_05.mat'))
end
disp('Dendrites state 05 superagglos loaded/generated')