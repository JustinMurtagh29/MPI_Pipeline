function createNewDendriteSuperagglos(param, state)
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
   
    dataDir = fullfile(param.saveFolder, 'aggloState');

    [~, suffix, ~, dendriteVersion, dendriteVersionNew] = ...
        connectEM.setDendriteQueryState(state);

    % Load current state of agglomerates
    agglos = load(fullfile( ...
        dataDir, strcat('dendrites_',dendriteVersion,'.mat')));
    superAgglos = agglos.dendrites(agglos.indBigDends);
    
    % NOTE(amotta): `agglos` contains all super-agglomerates. `superAgglos`
    % is restricted to the large-enough ones.

    % Load linkages
    flightPaths = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'startAgglo','endAgglo','ff');
    
    m = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'idxGood');
    idxGood = m.idxGood;
    clear m
    
    m = load(fullfile(dataDir, strcat('axonDendriteQueryOverlaps_',suffix,'.mat')),'idxNoClearEndAndNoAxonAttachment');
   
    % flight paths with no clear end and no axon attachment (dangling)
    % should be added
    idxGood(m.idxNoClearEndAndNoAxonAttachment) = true;
    clear m
    
    flightPaths.startAgglo = flightPaths.startAgglo(idxGood);
    flightPaths.endAgglo = flightPaths.endAgglo(idxGood);
    flightPaths.ff = structfun(@(x)x(idxGood), flightPaths.ff, 'uni', 0);

    linkagesAgglos = cellfun(@(x,y) connectEM.edgeCreator(x,y,1), flightPaths.startAgglo, flightPaths.endAgglo, 'uni', 0);
    linkagesLUT = repelem(1:numel(linkagesAgglos),cellfun(@(x) size(x,1),linkagesAgglos));
    linkagesAgglos = sort(cat(1, linkagesAgglos{:}),2);
    
    % reshape empty cells to 1st dimension
    flightPaths.endAgglo = cellfun(@(x) reshape(x,[],1),flightPaths.endAgglo,'uni',0);
    % sanity checks
%     assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
%     assert(size(linkagesAgglos, 1) == numel(cat(1,flightPaths.endAgglo{:})));
%     assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
    assert(~any(diff(structfun(@numel, flightPaths.ff))));
    
    % Eliminate duplicates
    [~, positionUniques] = unique(linkagesAgglos, 'rows');
    duplicates = true(size(linkagesAgglos,1),1);
    duplicates(positionUniques) = 0;
    linkagesAgglos = linkagesAgglos(~duplicates,:);
    
    % eliminate only those start/endAgglo/ffs where all found ends result
    % were found as duplicate
    duplicates = find(histc(linkagesLUT(duplicates),1:max(linkagesLUT)) == histc(linkagesLUT,1:max(linkagesLUT)));
    positionUniques = setdiff(1:max(linkagesLUT),duplicates);
    flightPaths.startAgglo(duplicates,:) = [];
    flightPaths.endAgglo(duplicates,:) = [];
    flightPaths.ff = structfun(@(x)x(positionUniques), flightPaths.ff, 'uni', 0);
    
    eqClassCCfull = Graph.findConnectedComponents(linkagesAgglos(~any(isnan(linkagesAgglos),2),:), true, true);
    eqClassCCfull = [eqClassCCfull; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCCfull)))'];
    
    dendritesNew = connectEM.mergeSuperagglosBasedOnFlightPath( ...
        superAgglos, eqClassCCfull, flightPaths.startAgglo, ...
        flightPaths.endAgglo, flightPaths.ff);

    % sanity check
    assert(all(arrayfun(@(x) ...
        numel(Graph.findConnectedComponents(x.edges)) == 1, ...
        dendritesNew(arrayfun(@(x) size(x.nodes, 1), dendritesNew) > 1))));

    eqClasses = flightPaths.eqClassCCfull;
    
    out = struct;
    out.param = param;
    out.state = state;
    out.gitInfo = Util.gitInfo();
    
    out.eqClasses = eqClasses;
    
    out.dendrites = cat(1, dendritesNew, agglos.dendrites(~agglos.indBigDends));
    out.indBigDends = false(length(out.dendrites),1);
    out.indBigDends(1:numel(dendritesNew), 1) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, sprintf('dendrites_%s.mat', dendriteVersionNew));
    Util.saveStruct(saveFile, out);
    system(sprintf('chmod -w "%s"', saveFile));
end
