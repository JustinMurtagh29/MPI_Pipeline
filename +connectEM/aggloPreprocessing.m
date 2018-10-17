%% load graph etc
% Comment MB: executed with clean state, e.g. ~exist will all trigger

overwrite = 0;
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');
info = Util.runInfo(); % added by BS
statesDendrites = {'dendrites_01','dendrites_02','dendrites_03_v2','dendrites_03_v2_splitmerged','dendrites_04','dendrites_05','dendrites_06','dendrites_07','dendrites_08','dendrites_09','dendrites_10','dendrites_11','dendrites_12','dendrites_13','dendrites_14','dendrites_15','dendrites_16','dendrites_16_b','dendrites_wholeCells_01_v5'};
statesAxons = {'axons_01','axons_02','axons_03'};
statesWC = {'wholeCells_01','wholeCells_02','wholeCells_03','wholeCells_04','wholeCells_05','wholeCells_06','wholeCells_07','wholeCells_08','wholeCells_GTAxon_08_v4'};
existentDendrites = cellfun(@(x) exist(fullfile(outputFolder,strcat(x,'.mat')),'file'),statesDendrites) | overwrite;
existentAxons = cellfun(@(x) exist(fullfile(outputFolder,strcat(x,'.mat')),'file'),statesAxons) | overwrite;
existentWC = cellfun(@(x) exist(fullfile(outputFolder,strcat(x,'.mat')),'file'),statesWC) | overwrite;

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
sMpoints(cat(1,segmentMeta.segIds),1:3) = cat(1,segmentMeta.point');
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
    %     [dendrites,dendriteLUT] = connectEM.applyWholeCellCorrections(dendrites,somaAgglos,p,fullfile(outputFolder,correctionFolder),1,sMpoints);
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),1,sMpoints);
    
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
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,sMpoints);
    
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
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_04');
    fprintf('Folder with correction nmls for dendrites_04 is %s\n',correctionFolder);
    %     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    
    % split the stuff to be added now
    dendrites = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    axons = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    
    save(fullfile(outputFolder,'wholeCells_01.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_05.mat'),'dendrites','WholeCellId','myelinDend','indBigDends','info');
elseif ~existentDendrites(7) || ~existentWC(2)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_05.mat'))
    load(fullfile(outputFolder,'wholeCells_01.mat'));
end
disp('Dendrites state 05 dendrites loaded/generated')

%% next round of manual whole cell fixes
if  ~existentWC(2) || ~existentDendrites(7)
    
    % this function is in Alessandro's repo
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    % ending detection done by christian!
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_05');
    fprintf('Folder with correction nmls for dendrites_05 is %s\n',correctionFolder);
    %     connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    % no correction was necessary
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    
    save(fullfile(outputFolder,'wholeCells_02.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_06.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(8) || ~existentWC(3)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_06.mat'))
    load(fullfile(outputFolder,'wholeCells_02.mat'),'wholeCells');
end
disp('State 06 dendrites loaded/generated')
%% next round of manual whole cell fixes
if ~existentWC(3) || ~existentDendrites(8)
    
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    % ending detection done by christian!
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_06');
    fprintf('Folder with correction nmls for dendrites_06 is %s\n',correctionFolder);
    %     connectEM.applyAggloCorrections(cat(1,wholeCells,dendrites),p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    % no correction was necessary
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    save(fullfile(outputFolder,'wholeCells_03.mat'),'wholeCells');
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

%% next round of manual whole cell fixes
if ~existentDendrites(9)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    [dendriteLUT,dendriteSegIds] = Superagglos.buildLUT(dendrites);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = unique(accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));
    
    % writes out skeletons of whole cells for inspection
    %     connectEM.aggloAutoView('dendrites_07','wc_border',1)
    
    correctionFolder = 'WholeCellCorrections_07';
    fprintf('Folder with correction nmls for state dendrites_07 is %s\n',fullfile(outputFolder,correctionFolder));
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),1,sMpoints);
    
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

%% next round of manual whole cell fixes
if ~existentDendrites(10)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls
    
    correctionFolder = 'WholeCellCorrections_08';
    fprintf('Folder with correction nmls for state dendrites_08 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
    %     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    
    % split the stuff to be added now
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,sMpoints,axons,1);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = unique(accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_09.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(11)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_09.mat'))
end
disp('State 09 dendrites loaded/generated')

%% next round of manual whole cell fixes
if ~existentDendrites(11)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls
    
    correctionFolder = 'WholeCellCorrections_09';
    fprintf('Folder with correction nmls for state dendrites_09 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
    %     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,sMpoints,axons,1);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = unique(accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_10.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(12)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_10.mat'))
end
disp('State 10 dendrites loaded/generated')


%% next round of manual whole cell fixes
if ~existentDendrites(12)
    load(fullfile(outputFolder,'axons_07_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls
    
    correctionFolder = 'WholeCellCorrections_10';
    fprintf('Folder with correction nmls for state dendrites_10 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
    %     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,sMpoints,axons,1);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = unique(accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_11.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(13)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_11.mat'))
end
disp('State 11 dendrites loaded/generated')

%% next round of manual whole cell fixes
if ~existentDendrites(13)
    
    load(fullfile(outputFolder,'axons_08_b.mat'),'axons')
    
    % here was christians ending detection done on dendrites_08 + manual
    % correction of nmls
    
    correctionFolder = 'WholeCellCorrections_11';
    fprintf('Folder with correction nmls for state dendrites_11 is %s\n',fullfile(outputFolder,correctionFolder));
    
    % write out the stuff to be added for checking
    %     connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),2,sMpoints,axons,1);
    
    % split the stuff to be added now (dendrites had no checks necessary)
    [axons,axonLUT] = connectEM.applyAggloCorrections(axons,p,fullfile(outputFolder,correctionFolder,'checkedBeforeAdd'),1,sMpoints);
    
    
    [dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendrites,p,fullfile(outputFolder,correctionFolder),0,sMpoints,axons,1);
    
    % test for duplets
    agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
    assert(numel(agglos)==numel(unique(agglos)))
    
    dendriteSegIds = find(dendriteLUT);
    [ismem,ind] = ismember(somaSegIds,dendriteSegIds);
    % get each dend id which contains most of the seg ids of each soma
    BorderWholeCellId = unique(accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));
    
    indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_12.mat'),'dendrites','BorderWholeCellId','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(14)
    clear dendrites myelinDend indBigDends
    load(fullfile(outputFolder,'dendrites_12.mat'))
end
disp('State 12 dendrites loaded/generated')

%% next round of manual whole cell fixes
if ~existentDendrites(14) || ~existentWC(4)
    load(fullfile(outputFolder,'axons_08_b.mat'),'axons')
    
    % the corrections are based on checking border whole cells in
    % dendrites_12 for remaining mergers/missing ends
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_13');
    fprintf('Folder with correction nmls for dendrites_12 is %s\n',correctionFolder);
    [ borderwholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    
    load(fullfile(outputFolder,'wholeCells_03.mat'),'wholeCells');
    wholeCells = cat(1,wholeCells,borderwholeCells);
    save(fullfile(outputFolder,'wholeCells_04.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_13.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(15) || ~existentWC(5)
    load(fullfile(outputFolder,'wholeCells_04.mat'))
    load(fullfile(outputFolder,'dendrites_13.mat'))
end
disp('State 13 dendrites loaded/generated')

%% load all soma whole cell agglos
[somaAgglos,somaCoords,KAMINcoords] = connectEM.getSomaAgglos(fullfile(outputFolder,'somas_with_merged_somas.mat'),'all');
somaCoords = somaCoords(~cellfun(@isempty,somaAgglos)); % remove not found somata
somaAgglos = somaAgglos(~cellfun(@isempty,somaAgglos)); % remove not found somata
somaSegIds = cell2mat(somaAgglos);
% remove duplicate segIds
[~,ic] = unique(somaSegIds);
duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
% somaDuplicateIds = cellfun(@(x) any(intersect(x,duplicates)),somaAgglos);
[somaAgglos,ia] = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
somaCoords = arrayfun(@(x) somaCoords{x}(ia{x},:),1:numel(ia),'uni',0);
somaSegIds = cell2mat(somaAgglos);
somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));
disp('Soma whole cell agglos loaded')

%% next round of manual whole cell fixes
if  ~existentWC(5) || ~existentDendrites(15)
    % this script was used to write out whole cells that mighted have
    % skipped agglos inbetween
    %     connectEM.searchSkippedAgglos(wholeCells,'/tmpscratch/mbeining/checkForMissingEndings',p)
    
    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    
    % the corrections are based on checking the output of
    % connectEM.searchSkippedAgglos
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_14');
    fprintf('Folder with correction nmls for dendrites_13 is %s\n',correctionFolder);
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    
    save(fullfile(outputFolder,'wholeCells_05.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_14.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(16) || ~existentWC(6)
    load(fullfile(outputFolder,'wholeCells_05.mat'))
    load(fullfile(outputFolder,'dendrites_14.mat'))
end
disp('State 14 dendrites loaded/generated')

%% next round of manual whole cell fixes
if  ~existentWC(6) || ~existentDendrites(16)
    % created nmls with connectEM.autoView('wholeCells_05','cells') and
    % checked only the five corrected last time
    
    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    
    % the corrections are based on checking the output of
    % connectEM.searchSkippedAgglos
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_15');
    fprintf('Folder with correction nmls for dendrites_14 is %s\n',correctionFolder);
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    
    save(fullfile(outputFolder,'wholeCells_06.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_15.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(17) || ~existentWC(7)
    load(fullfile(outputFolder,'wholeCells_06.mat'))
    load(fullfile(outputFolder,'dendrites_15.mat'))
end
disp('State 15 dendrites loaded/generated')

%% next round of manual whole cell fixes
if  ~existentWC(7) || ~existentDendrites(17)
    % created nmls with connectEM.autoView('wholeCells_05','cells') and
    % checked only the five corrected last time
    
    load(fullfile(outputFolder,'axons_09_a.mat'),'axons')
    correctionFolder = fullfile(outputFolder,'WholeCellCorrections_16');
    fprintf('Folder with correction nmls for dendrites_15 is %s\n',correctionFolder);
    [ wholeCells,dendrites,indBigDends,myelinDend ] = connectEM.wcApplyManualAnnotation(cat(1,wholeCells,dendrites), correctionFolder, axons,somaSegIds,somaLUT, sMpoints );
    save(fullfile(outputFolder,'wholeCells_07.mat'),'wholeCells');
    save(fullfile(outputFolder,'dendrites_16.mat'),'dendrites','myelinDend','indBigDends')%,'info');
elseif ~existentDendrites(18) || ~existentWC(8) || ~existentWC(9)
    load(fullfile(outputFolder,'wholeCells_07.mat'))
    load(fullfile(outputFolder,'dendrites_16.mat'))
end
disp('State 16 dendrites loaded/generated')

%% patch in dendrite dendrite and dangling flight paths
%{
if ~existentDendrites(18)
    connectEM.getDendriteQueryOverlapB(p,'2.3')
    connectEM.getDendQueryAxonAggloOverlapB(p,'2.3')
    connectEM.createNewDendriteSuperagglos(p,'2.3')
    load(fullfile(outputFolder,'dendrites_16_b.mat'))
end
%% patch in AIS

%}
%% remove axons from whole cells and add them to a new state of dendrites

% assure that wholeCells are sorted in the same way as the somaAgglos
newInd = arrayfun(@(x) mode(nonzeros(somaLUT(x.nodes(~isnan(x.nodes(:,4)),4)))),wholeCells);
assert(numel(unique(newInd)) == numel(newInd))
[~,newInd] = sort(newInd);
wholeCells = wholeCells(newInd);

if  ~existentWC(8)
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
        %     branchLUT = Superagglos.buildLUT(branches);
        %     countPre = histc(branchLUT(presynSegIds),0:numel(branches));
        %     countPost = histc(branchLUT(postsynSegIds),0:numel(branches));
        
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
        %         tmp = rmfield(Superagglos.removeNodesFromAgglo(wholeCells(n),~wholeCells(n).axon),'axon'); % get the whole cells without the axon
        %         ix = max(arrayfun(@(x) size(x.nodes,1),tmp));
        %         tmp = tmp(ix);
        %         AIS(n) = tmp;
    end
    % remove all branches that have been labeled axonic from all wholeCells
    
    % concatenate truncated whole cells with dendrite class and make new state
    indWholeCells = cat(1,false(numel(dendrites),1),true(numel(wholeCellsNoAxon),1));
    dendrites = cat(1,dendrites,wholeCellsNoAxon');
    indBigDends = cat(1,indBigDends,true(numel(wholeCellsNoAxon),1));
    indAIS = cat(1,indAIS,false(numel(wholeCellsNoAxon),1));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    
    save(fullfile(outputFolder,'wholeCells_08.mat'),'wholeCells');
    
    save(fullfile(outputFolder,'dendrites_wholeCells_01.mat'),'dendrites','myelinDend','indBigDends','indWholeCells','indAIS')%,'info');
end

if  ~existentWC(9)
    % now comes the manual axon detection part
    clear wholeCellsNoAxon
    
    dirWCGT = {'/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/centerWholeCellGT/axonInfo','/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/borderWholeCellGT_axonId'};
    filenames = cell(0);
    for a = 1:numel(dirWCGT)
        files = dir(fullfile(dirWCGT{a},'*.nml'));
        filenames = cat(2,filenames,cellfun(@(x) fullfile(dirWCGT{a},x),{files.name},'uni',0));
    end
    for n = 1:numel(wholeCells)  % make the axon field NaN for all wholeCells, to have those wCs labeled which do not have an axon marked in the GT
        wholeCells(n).axon = NaN(size(wholeCells(n).nodes,1),1);
    end
    wcLUT = Superagglos.buildLUT(wholeCells,segmentMeta.maxSegId);
    
    usedCells = zeros(numel(wholeCells),1);
    skelsNotFound = [];
    for f = 1:numel(filenames)
        skel = skeleton(filenames{f});
        PL(f) = skel.pathLength/1000;
        numSkelNodes = cellfun(@(x) size(x,1),skel.nodes);
        skelLUT = repelem(1:numel(skel.nodes),numSkelNodes);
        warning('off')
        skelSegIds = Seg.Global.getSegIds(p,cell2mat(cellfun(@(x) x(:,1:3),skel.nodes,'uni',0)));  % extract the seg Ids of all skel nodes
        warning on
        ind = mode(nonzeros(wcLUT(nonzeros(skelSegIds)))); % get the whole cell overlapping the most with the skeleton in terms of segIds
        if isnan(ind) || sum(wcLUT(nonzeros(skelSegIds))==ind)/numel(nonzeros(skelSegIds)) < 0.33 % if no overlap found or overlap is less than one third of the skeleton
            warning('Found no corresponding whole Cell to skeleton from file %s. Trying to use somaAgglo as index...',filenames{f})
            ind = mode(nonzeros(somaLUT(nonzeros(skelSegIds))));
            if isnan(ind)
                % last try to retrieve correct whole cell by checking
                % neighbours of the skel, but only if skel is very small
                % (as it would be for a border soma)
                if sum(numSkelNodes) > 10 || isnan(mode(nonzeros(somaLUT(nonzeros(cat(1,graph.neighbours{nonzeros(skelSegIds)}))))))
                    warning('Still no corresponding whole Cell to skeleton from file %s found! Skipping this one...',filenames{f})
                    skelsNotFound = cat(1,skelsNotFound,f);
                    continue
                else
                    ind = mode(nonzeros(somaLUT(nonzeros(cat(1,graph.neighbours{nonzeros(skelSegIds)})))));
                end
            end
        end
        if usedCells(ind)
            error('Error skeleton %s: Cell %d already used from skeleton %s!',filenames{f},ind,filenames{usedCells(ind)})
        else
            usedCells(ind) = f;
        end
        
        if numel(numSkelNodes) > 1  % if there is only one skeleton, there was no axon found
            indAxonSkel = find(~cellfun(@isempty,regexp(skel.names,'axon')));
            if isempty(indAxonSkel)  % if skeleton was not marked with an axon label, check comments for axon label
                [comments,treeIndComment] = skel.getAllComments;
                indComment = ~cellfun(@isempty,regexp(comments,'axon','ONCE'));
                indAxonSkel = unique(treeIndComment(indComment));
                if isempty(indAxonSkel)
                    % if also nothing found in comments, take the smaller one of the two skeletons as axon
                    [~,indAxonSkel] = min(numSkelNodes);
                end
            end
            % get node indices in whole cell which overlap with GT
            [~,mask] = ismember(unique(nonzeros(skelSegIds(skelLUT==indAxonSkel))),wholeCells(ind).nodes(:,4));
            mask = sort(unique(nonzeros(mask)));
            centerSoma = mean(somaCoords{ind},1); % get center coordinate of soma
            % get distance of all nodes to center
            distToCenter = sqrt(sum(bsxfun(@minus,bsxfun(@times,wholeCells(ind).nodes(:,1:3),[11.24,11.24,28]),centerSoma).^2,2));
            % get the minimal distance of the axon branch to the soma plus
            % some threshold
            minDist = min(distToCenter(mask))+2000;
            if minDist < 2000
                warning('Distance of axon to center of soma is less than 2 micron. Seems strange.. (%s)',filenames{f})
            end
            % grow out soma branch to get all nodes of the branch but do
            % not go below the minDist threshold
            while true
                maskNew = unique(wholeCells(ind).edges(any(ismember(wholeCells(ind).edges,mask),2),:));
                if isequal(mask,maskNew(distToCenter(maskNew) > minDist))
                    break
                else
                    mask = maskNew(abs(distToCenter(maskNew)) > minDist);
                end
            end
            wholeCells(ind).axon = false(size(wholeCells(ind).nodes,1),1);
            wholeCells(ind).axon(mask) = true;
        end
    end
    fprintf('Total path length of all whole Cells (GT based) is %g micron\n',sum(PL(usedCells)));
    if ~isempty(skelsNotFound)
        warning('Skeletons of files \n%s\n did not have a whole cell partner!',filenames{skelsNotFound})
    end
    undefAxon = arrayfun(@(x) any(isnan(x.axon)),wholeCells);
    nodesToDelete = arrayfun(@(x) find(x.axon),wholeCells,'uni',0);
    nodesToDelete(undefAxon) = cell(sum(undefAxon),1);
    % stupid splits can occur during axon deletion when a branch node was
    % the only connection between a segID being in somaAgglo and one (or more)
    % being not, so check each deletion and choose the biggest agglo as the
    % one to take as wholeCell...
    clear wholeCellsNoAxon
    for n = 1:numel(wholeCells)
        tmp = Superagglos.removeNodesFromAgglo(wholeCells(n),nodesToDelete(n));
        if numel(tmp)>1
            [~,ind] = max(arrayfun(@(x) size(x.nodes,1),tmp));
        else
            ind = 1;
        end
        wholeCellsNoAxon(n) = rmfield(tmp(ind),'axon');
    end
    save(fullfile(outputFolder,'wholeCells_GTAxon_08_v4.mat'),'wholeCells');
    
    load(fullfile(outputFolder,'dendrites_16.mat'))
    % concatenate truncated whole cells with dendrite class and make new state
    indWholeCells = cat(1,false(numel(dendrites),1),true(numel(wholeCellsNoAxon),1));
    dendrites = cat(1,dendrites,wholeCellsNoAxon');
    indBigDends = cat(1,indBigDends,true(numel(wholeCellsNoAxon),1));
    indAIS = cat(1,indAIS,false(numel(wholeCellsNoAxon),1));
    [ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    save(fullfile(outputFolder,'dendrites_wholeCells_01_v5.mat'),'dendrites','myelinDend','indBigDends','indWholeCells','indAIS')%,'info');
    
end
%%
