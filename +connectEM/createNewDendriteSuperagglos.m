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
    flightPaths = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'edgesCC','eqClassCCfull','startAgglo','endAgglo','ff');
    
    m = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'idxGood');
    idxGood = m.idxGood;
    clear m
    
    m = load(fullfile(dataDir, strcat('axonDendriteQueryOverlaps_',suffix,'.mat')),'idxNoClearEndButAxonAttachment');
    idxNoClearEndButAxonAttachment = m.idxNoClearEndButAxonAttachment;
    clear m
    % remove flight paths with no clear and but clear axon attachment
    idxGood(idxNoClearEndButAxonAttachment) = false;

    flightPaths.startAgglo = flightPaths.startAgglo(idxGood);
    flightPaths.endAgglo = flightPaths.endAgglo(idxGood);
    flightPaths.ff = structfun(@(x)x(idxGood), flightPaths.ff, 'uni', 0);

    linkagesAgglos = cellfun(@(x,y) connectEM.edgeCreator(x,y,1), flightPaths.startAgglo, flightPaths.endAgglo, 'uni', 0);
    linkagesLUT = repelem(1:numel(linkagesAgglos),cellfun(@(x) size(x,1),linkagesAgglos));
    linkagesAgglos = cat(1, linkagesAgglos{:});
    
    % reshape empty cells to 1st dimension
    flightPaths.endAgglo = cellfun(@(x) reshape(x,[],1),flightPaths.endAgglo,'uni',0);
    % sanity checks
%     assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
%     assert(size(linkagesAgglos, 1) == numel(cat(1,flightPaths.endAgglo{:})));
%     assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
    assert(~any(diff(structfun(@numel, flightPaths.ff))));
    
    % Eliminate duplicates
    [~, positionUniques] = unique(sort(linkagesAgglos,2), 'rows');
    duplicates = true(size(linkagesAgglos,1),1);
    duplicates(positionUniques) = 0;
    
    % NOTE(amotta): Dimensions agree!
    flightPaths.startAgglo(linkagesLUT(duplicates),:) = [];
    flightPaths.endAgglo(linkagesLUT(duplicates),:) = [];
    flightPaths.ff = structfun(@(x)x(linkagesLUT(~duplicates)), flightPaths.ff, 'uni', 0);
    
    dendritesNew = connectEM.mergeSuperagglosBasedOnFlightPath( ...
        superAgglos, flightPaths.eqClassCCfull, flightPaths.startAgglo, ...
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
