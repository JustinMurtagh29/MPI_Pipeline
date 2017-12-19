%% load graph etc
% Comment MB: executed with clean state, e.g. ~exist will all trigger
overwrite = 0;
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');
info = Util.runInfo(); % added by BS
statesDendrites = {'dendrites_01','dendrites_02','dendrites_03_v2','dendrites_03_v2_splitmerged','dendrites_04','dendrites_05','dendrites_06','dendrites_07','dendrites_08','dendrites_09','dendrites_10','dendrites_11','dendrites_12','dendrites_13','dendrites_14','dendrites_15','dendrites_16','dendrites_andWholeCells_01'};
statesAxons = {'axons_01','axons_02','axons_03'};
statesWC = {'wholeCells_01','wholeCells_02','wholeCells_03','wholeCells_04','wholeCells_05','wholeCells_06','wholeCells_07'};
existentDendrites = cellfun(@(x) exist(strcat(x,'.mat'),'file'),statesDendrites) | overwrite;
existentAxons = cellfun(@(x) exist(strcat(x,'.mat'),'file'),statesAxons) | overwrite;
existentWC = cellfun(@(x) exist(strcat(x,'.mat'),'file'),statesWC) | overwrite;

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
if ~existentDendrites(1)
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
    save(fullfile(outputFolder,'dendritesEdgesGTall.mat'), 'edgesGTall','info');
    
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
    save(fullfile(outputFolder,'dendrites_01.mat'),'dendrites','indSingleSeg','info');
elseif ~existentDendrites(2)
    load(fullfile(outputFolder,'dendrites_01.mat'),'dendrites');
end

%% load/ create axon superagglos
if ~existentAxons(1)
    % all axons agglomerate results from after directionality based grid
    % search (> 100 voxel)
    load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');
    load('/tmpscratch/mberning/edgesGTall.mat', 'edgesGTall');
    % remove all duplicates from hot edge list
    edgesGTall = unique(sort(edgesGTall,2),'rows');
    save(fullfile(outputFolder,'axonsEdgesGTall.mat'), 'edgesGTall','info');
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
    save(fullfile(outputFolder,'axons_01.mat'),'axons','info');
elseif ~existentAxons(2)
    load(fullfile(outputFolder,'axons_01.mat'),'axons');
end
disp('State 01 superagglos loaded/generated')

%% execute corresponding edges
if ~existentAxons(2)
    axons = connectEM.executeEdges(axons,corrEdges,segmentMeta);
    save(fullfile(outputFolder,'axons_02.mat'),'axons','info');
elseif ~existentAxons(3)
    load(fullfile(outputFolder,'axons_02.mat'),'axons');
end
if ~existentDendrites(2)
    dendrites = connectEM.executeEdges(dendrites,corrEdges,segmentMeta);
    indSingleSeg = arrayfun(@(x) size(x.nodes,1),dendrites)==1;
    save(fullfile(outputFolder,'dendrites_02.mat'),'dendrites','indSingleSeg','info');
elseif ~existentDendrites(3)
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
if ~existentAxons(3) || ~existentDendrites(3)
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
    
    save(fullfile(outputFolder,'axons_03.mat'),'axons','myelinAxon','indBigAxons','info');
    save(fullfile(outputFolder,'dendrites_03_v2.mat'),'dendrites','myelinDend','indBigDends','info');
elseif ~existentDendrites(4)
%     load(fullfile(outputFolder,'axons_03.mat'),'axons','myelinAxon','indBigAxons');
    load(fullfile(outputFolder,'dendrites_03_v2.mat'),'dendrites','myelinDend','indBigDends');
end
disp('State 03 superagglos loaded/generated')

%% load soma whole cell agglos
somaAgglos = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),'center');
somaAgglos = somaAgglos(~cellfun(@isempty,somaAgglos)); % remove not found somata
somaSegIds = cell2mat(somaAgglos);
% remove duplicate segIds
[~,ic] = unique(somaSegIds);
duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
somaAgglos = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
somaSegIds = cell2mat(somaAgglos);
somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
disp('Soma whole cell agglos loaded')


%% apply manual fixes of whole cells in dendrite class
if ~existentDendrites(4)
    % ending detection done by christian!
    correctionFolder = 'WholeCellCorrections_03_v2';
    fprintf('Folder with correction nmls for state dendrites_03_v2 is %s\n',fullfile(outputFolder,correctionFolder));
%     [dendrites,dendriteLUT] = connectEM.applyWholeCellCorrections(dendrites,somaAgglos,p,fullfile(outputFolder,correctionFolder),1);
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),1);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_03_v2_splitmerged.mat'),'dendrites','WholeCellId','myelinDend','indBigDends','info');
elseif ~existentDendrites(5)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_03_v2_splitmerged.mat'),'dendrites','WholeCellId')%,'myelinDend','indBigDends');
end
disp('Dendrites state 03 splitmerged superagglos loaded/generated')

%% next round of fixes
% apply manual fixes of whole cells in dendrite class
if ~existentDendrites(5)
    % ending detection done by christian!
    correctionFolder = 'WholeCellCorrections_03_v2_splitmerged';
    fprintf('Folder with correction nmls for state dendrites_03_v2_splitmerged is %s\n',fullfile(outputFolder,correctionFolder));
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder));

    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_04.mat'),'dendrites','WholeCellId','myelinDend','indBigDends','info');
elseif ~existentDendrites(6)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_04.mat'))
end
disp('Dendrites state 04 dendrites loaded/generated')
%% next round of fixes
% apply manual fixes of whole cells in dendrite class
if ~existentDendrites(6) || ~existentWC(1)
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    % ending detection done by christian!
    correctionFolder = 'WholeCellCorrections_04';
    fprintf('Folder with correction nmls for state dendrites_04 is %s\n',fullfile(outputFolder,correctionFolder));
%     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now
    dendrites = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    wholeCells = dendrites(WholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),WholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_01.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_05.mat'),'dendrites','WholeCellId','myelinDend','indBigDends','info');
    
elseif ~existentDendrites(7) || ~existentWC(2)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_05.mat'))
    load(fullfile(outputFolder,'wholeCells_01.mat'));
end
disp('Dendrites state 05 dendrites loaded/generated')

%% 
if  ~existentWC(2) || ~existentDendrites(7)
    
    % this function is in Alessandro's repo
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    % ending detection done by christian!
    correctionFolder = 'WholeCellCorrections_05';
    fprintf('Folder with correction nmls for state dendrites_05 is %s\n',fullfile(outputFolder,correctionFolder));
%     connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),2,axons);
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    wholeCells = dendrites(WholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),WholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_02.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_06.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(8) || ~existentWC(3)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_06.mat'))
    load(fullfile(outputFolder,'wholeCells_02.mat'),'wholeCells');
end
disp('State 06 dendrites loaded/generated')
%%
if ~existentWC(3) || ~existentDendrites(8)
    
    % this function is in Alessandro's repo
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    % ending detection done by christian!
    correctionFolder = 'WholeCellCorrections_06';
    fprintf('Folder with correction nmls for state dendrites_06 is %s\n',fullfile(outputFolder,correctionFolder));
%     connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),2,axons);
    % no correction was necessary
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    WholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);
 
    wholeCells = dendrites(WholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),WholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_03.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_07.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(9) || ~existentWC(4)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_07.mat'))
    load(fullfile(outputFolder,'wholeCells_03.mat'),'wholeCells');
end
disp('State 07 dendrites loaded/generated')

%% load border soma whole cell agglos
somaAgglos = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),'border');
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

%% 
if ~existentDendrites(9)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    [dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);

    % writes out skeletons of whole cells for inspection
%     connectEM.aggloAutoView('dendrites_07','wc_border',1)
    
    correctionFolder = 'WholeCellCorrections_07';
    fprintf('Folder with correction nmls for state dendrites_07 is %s\n',fullfile(outputFolder,correctionFolder));
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),1);
  
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_08.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(10)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_08.mat'))
end
disp('State 08 dendrites loaded/generated')

%% 
if ~existentDendrites(10)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls   
    
    correctionFolder = 'WholeCellCorrections_08';
    fprintf('Folder with correction nmls for state dendrites_08 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
%     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now
    dendrites = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);

    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       
  
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_09.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(11)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_09.mat'))
end
disp('State 09 dendrites loaded/generated')

%% 
if ~existentDendrites(11)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls   
    
    correctionFolder = 'WholeCellCorrections_09';
    fprintf('Folder with correction nmls for state dendrites_09 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
%     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);

    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       
  
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_10.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(12)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_10.mat'))
end
disp('State 10 dendrites loaded/generated')


%% 
if ~existentDendrites(12)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls   
    
    correctionFolder = 'WholeCellCorrections_10';
    fprintf('Folder with correction nmls for state dendrites_10 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
%     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       
  
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_11.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(13)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_11.mat'))
end
disp('State 11 dendrites loaded/generated')

%% 
if ~existentDendrites(13)
    
    load(fullfile(outputFolder,'axons_08_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls   
    
    correctionFolder = 'WholeCellCorrections_11';
    fprintf('Folder with correction nmls for state dendrites_11 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
%     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,axons);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       
  
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_12.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(14)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_12.mat'))
end
disp('State 12 dendrites loaded/generated')

%% 
if ~existentDendrites(14) || ~existentWC(4)
    load(fullfile(outputFolder,'axons_08_b.mat'),'axons')
    
    % the corrections are based on checking border whole cells in
    % dendrites_12 for remaining mergers/missing ends
    correctionFolder = 'WholeCellCorrections_13';
    fprintf('Folder with correction nmls for state dendrites_12 is %s\n',fullfile(outputFolder,correctionFolder));
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,axons);

    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       

    load(fullfile(outputFolder,'wholeCells_03.mat'),'wholeCells');
    wholeCells = cat(1,wholeCells,dendrites(BorderWholeCellId));
    dendrites = dendrites(setdiff(1:numel(dendrites),BorderWholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_04.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_13.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(15) || ~existentWC(5)
    load(fullfile(outputFolder,'wholeCells_04.mat'))
    load(fullfile(outputFolder,'dendrites_13.mat'))
end
disp('State 13 dendrites loaded/generated')

%% load all soma whole cell agglos
somaAgglos = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),'all');
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

%% check manually added edges if agglos have been skipped there

if  ~existentWC(5) || ~existentDendrites(15)
    % this script was used to write out whole cells that mighted have
    % skipped agglos inbetween
%     connectEM.searchSkippedAgglos(wholeCells,'/tmpscratch/mbeining/checkForMissingEndings',p)

    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    
    % the corrections are based on checking the output of
    % connectEM.searchSkippedAgglos
    correctionFolder = 'WholeCellCorrections_14';
    fprintf('Folder with correction nmls for state dendrites_13 is %s\n',fullfile(outputFolder,correctionFolder));
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    wholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       

    wholeCells = dendrites(wholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),wholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_05.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_14.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(16) || ~existentWC(6)
    load(fullfile(outputFolder,'wholeCells_05.mat'))
    load(fullfile(outputFolder,'dendrites_14.mat'))
end
disp('State 14 dendrites loaded/generated')

%% another round of fixing added stuff

if  ~existentWC(6) || ~existentDendrites(16)
    % created nmls with connectEM.autoView('wholeCells_05','cells') and
    % checked only the five corrected last time

    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    
    % the corrections are based on checking the output of
    % connectEM.searchSkippedAgglos
    correctionFolder = 'WholeCellCorrections_15';
    fprintf('Folder with correction nmls for state dendrites_14 is %s\n',fullfile(outputFolder,correctionFolder));
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    wholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       

    wholeCells = dendrites(wholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),wholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_06.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_15.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(17) || ~existentWC(7)
    load(fullfile(outputFolder,'wholeCells_06.mat'))
    load(fullfile(outputFolder,'dendrites_15.mat'))
end
disp('State 15 dendrites loaded/generated')

%% another round of fixing added stuff

if  ~existentWC(7) || ~existentDendrites(17)
    % created nmls with connectEM.autoView('wholeCells_05','cells') and
    % checked only the five corrected last time

    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    
    % the corrections are based on checking the output of
    % connectEM.searchSkippedAgglos
    correctionFolder = 'WholeCellCorrections_16';
    fprintf('Folder with correction nmls for state dendrites_15 is %s\n',fullfile(outputFolder,correctionFolder));
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),0,axons);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    wholeCellId = accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode);       

    wholeCells = dendrites(wholeCellId);
    dendrites = dendrites(setdiff(1:numel(dendrites),wholeCellId));
    
    save(fullfile(outputFolder,'wholeCells_07.mat'),'wholeCells');
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_16.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(18)
    load(fullfile(outputFolder,'wholeCells_07.mat'))
    load(fullfile(outputFolder,'dendrites_16.mat'))
end
disp('State 16 dendrites loaded/generated')
% %% add spines to all agglos
% if ~exist(fullfile(outputFolder,'wholeCells_04.mat'),'file') || overwrite
%     % this function is in Alessandro's repo
%     wholeCells = L4.Spine.Head.attachRun('ex145_07x2_roi2017','wholeCells_03.mat',1);
%     save(fullfile(outputFolder,'wholeCells_04.mat'),'wholeCells') % dendrites should not become much bigger or get more myelin when spines are attached
% else
%     load(fullfile(outputFolder,'wholeCells_04.mat'))
% end
% 
% 
% if ~exist(fullfile(outputFolder,'dendrites_08.mat'),'file') || overwrite
%     % this function is in Alessandro's repo
%     dendrites = L4.Spine.Head.attachRun('ex145_07x2_roi2017','dendrites_07.mat',1);
%     save(fullfile(outputFolder,'dendrites_08.mat'),'dendrites','myelinDend','indBigDends') % dendrites should not become much bigger or get more myelin when spines are attached
% else
%     clear dendrites myelinDend indBigDends
%     load(fullfile(outputFolder,'dendrites_08.mat'))
% end


%% remove axons from whole cells and add them to a new state of dendrites
assert(isequal(uint32(1:max(segmentMeta.segIds))',segmentMeta.segIds))
% load(fullfile(p.saveFolder,'connectomeState','SynapseAgglos_v2.mat'),'synapses') % (erzeugt via E:\workspace\pipeline\Benedikt\+L4\+Synapses\+Scripts\synapseDetection.m)
% presynSegIds = cellfun(@(x) x(1),synapses.presynId); % only one presyn id necessary
% postsynSegIds = cellfun(@(x) x(1),synapses.postsynId); % only one postsyn id necessary
% maxSegId = max([presynSegIds;postsynSegIds]);
clear wholeCellsNoAxon
for n = 1:numel(wholeCells)
    % get cell branches by removing soma agglo and small segments
    branches = Superagglos.removeSegIdsFromAgglos(wholeCells(n),somaAgglos{n});
    branches = branches(arrayfun(@(x) size(x.nodes,1),branches)>20); % remove all small leftovers smaller than 20 nodes
    branchLUT = Superagglos.buildLUT(branches,maxSegId);
    countPre = histc(branchLUT(presynSegIds),0:numel(branches));
    countPost = histc(branchLUT(postsynSegIds),0:numel(branches));

    if isempty(branches) % if there was only somatic stuff, nothing is axonic
        wholeCells(n).axon = false(size(wholeCells(n).nodes,1),1);
    else
        % get the median axon and dendrite probabilities for each branch
        [axonP,dendriteP] = arrayfun(@(x) deal(median(segmentMeta.axonProb(x.nodes(~isnan(x.nodes(:,4)),4))),median(segmentMeta.dendriteProb(x.nodes(~isnan(x.nodes(:,4)),4)))),branches);
%         [spineP,gliaP] = arrayfun(@(x) deal(median(segmentMeta.spineProb(x.nodes(~isnan(x.nodes(:,4)),4))),median(segmentMeta.gliaProb(x.nodes(~isnan(x.nodes(:,4)),4)))),branches);
        [~,ind1] = min(dendriteP);
        [~,ind2] = max(axonP);
        
        if ~any(axonP > 0.1) % one axon probability has at least be over 10% otherwise do not label any branch as axonic
            wholeCells(n).axon = false(size(wholeCells(n).nodes,1),1);
        elseif ind1==ind2 %if one branch has the highest axon and the lowest dendrite prob, then everything is easy
            wholeCells(n).axon = ismember(wholeCells(n).nodes(:,1:3),branches(ind1).nodes(:,1:3),'rows');
        else
            % check which ratio (dendriteProb or axonProb versus median) is
            % more prominent and accept only if axonProb of one branch is
            % high
            [~,ind] = max([median(dendriteP)/min(dendriteP), max(axonP)/median(axonP)]);
            switch ind
                case 1
                    error('No branch showed a high axon probability compared to the median of all branches.')
                case 2
                    % label the found branch as axonic
                    wholeCells(n).axon = ismember(wholeCells(n).nodes(:,1:3),branches(ind2).nodes(:,1:3),'rows');
            end
        end
    end
    tmp = rmfield(Superagglos.removeNodesFromAgglo(wholeCells(n),wholeCells(n).axon),'axon'); % get the whole cells without the axon
    tmp = tmp(arrayfun(@(x) size(x.nodes,1),tmp)>10);
    if numel(tmp) > 1
        [~,ind] = max(arrayfun(@(x) sum(ismember(somaAgglos{n},x.nodes(:,4))),tmp));
        wholeCellsNoAxon(n) = tmp(ind);
        warning('Multiple CCs found when deleting axonic branch of whole Cell %d. Please check',n)
    else
        wholeCellsNoAxon(n) = tmp;
    end
end
% remove all branches that have been labeled axonic from all wholeCells

% concatenate truncated whole cells with dendrite class and make new state
indWholeCells = cat(1,false(numel(dendrites),1),true(numel(wholeCellsNoAxon),1));
dendrites = cat(1,dendrites,wholeCellsNoAxon');
indBigDends = cat(1,indBigDends,true(numel(wholeCellsNoAxon),1));
[ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class

save(fullfile(outputFolder,'dendrites_andWholeCells_01.mat'),'dendrites','myelinDend','indBigDends','indWholeCells')%,'info');
%%
connectEM.getDendriteQueryOverlapB(p,'2.2')
connectEM.getDendQueryAxonAggloOverlapB(p,'2.2')
connectEM.createNewDendriteSuperagglos(p,'2.2')
