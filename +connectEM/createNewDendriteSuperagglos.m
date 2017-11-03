function createNewDendriteSuperagglos(param, state)
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
   
    dataDir = fullfile(param.saveFolder, 'aggloState');

    [~, suffix, dendriteVersion, dendriteVersionNew] = ...
        connectEM.setDendriteQueryState(state);

    % Load current state of agglomerates
    agglos = load(fullfile( ...
        dataDir, strcat('dendrites_',dendriteVersion,'.mat')));
    superAgglos = agglos.dendrites(agglos.indBigDends);
    
    % NOTE(amotta): `agglos` contains all super-agglomerates. `superAgglos`
    % is restricted to the large-enough ones.

    % Load linkages
    flightPaths = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'edgesCC','eqClassCCfull','startAgglo','endAgglo','ff');
    
    m = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')),'idxNoClearEnd','idxNoClearStart');
    idxNoClearEnd = m.idxNoClearEnd;
    idxNoClearStart = m.idxNoClearStart;
    idxNoClearStart = find(idxNoClearStart);
    idxNoClearEnd = find(idxNoClearEnd);
    
    clear m
    
    m = load(fullfile(dataDir, strcat('axonDendriteQueryOverlaps_',suffix,'.mat')),'idxNoClearAxonAttachment');
    idxNoClearAxonAttachment = m.idxNoClearAxonAttachment;
    clear m

    excludedFlights = idxNoClearEnd(~idxNoClearAxonAttachment);
    excludedFlights = cat(1,excludedFlights,idxNoClearStart);
%     endings = load(fullfile(dataDir, 'dendriteEndings.mat'));
%     endingOverlap = load(fullfile(dataDir, strcat('dendriteEndingOverlaps',suffix,'.mat')));
%     [linkages, idxMultipleHits] = cellfun(@connectEM.edgeCreator, endingOverlap.startEndingOverlaps, endingOverlap.endEndingOverlaps, 'uni', 0);
%     linkagesFlat = cat(1, linkages{:});
%     idxMultipleHits = cat(1, idxMultipleHits{:});
    
    linkagesAgglos = cellfun(@connectEM.edgeCreator, flightPaths.startAgglo, flightPaths.endAgglo, 'uni', 0);
    linkagesAgglos = cat(1, linkagesAgglos{:});
    
    idxExcludedFlights = true(size(flightPaths.startAgglo));
    idxExcludedFlights(excludedFlights) = false;
    flightPaths.startAgglo(excludedFlights,:) = [];
    flightPaths.endAgglo(excludedFlights,:) = [];
    flightPaths.ff = structfun(@(x)x(idxExcludedFlights), flightPaths.ff, 'uni', 0);
    linkagesAgglos(excludedFlights, :) = [];

    % sanity checks
    assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.endAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
    assert(~any(diff(structfun(@numel, flightPaths.ff))));
    
    % Eliminate duplicates
%     [~, positionUniques] = unique(sort(linkagesAgglos,2), 'rows');
%     duplicates = true(size(linkagesAgglos,1),1);
%     duplicates(positionUniques) = 0;
%     
%     % NOTE(amotta): Dimensions agree!
%     flightPaths.startAgglo(duplicates,:) = [];
%     flightPaths.endAgglo(duplicates,:) = [];
%     flightPaths.ff = structfun(@(x)x(~duplicates), flightPaths.ff, 'uni', 0);
    
    dendritesNew = connectEM.mergeSuperagglosBasedOnFlightPath( ...
        superAgglos, flightPaths.eqClassCCfull, flightPaths.startAgglo, ...
        flightPaths.endAgglo, flightPaths.ff);

    % sanity check
    assert(all(arrayfun(@(x) ...
        numel(Graph.findConnectedComponents(x.edges)) == 1, ...
        dendritesNew(arrayfun(@(x) size(x.nodes, 1), dendritesNew) > 1))));

    out = struct;
    out.param = param;
    out.state = state;
    out.gitInfo = Util.gitInfo();
    
    out.eqClasses = eqClasses;
    out.flightGroups = flightGroups;
    
    out.axons = cat(1, dendritesNew, agglos.dendrites(~agglos.indBigDends));
    out.indBigDends = false(length(out.axons),1);
    out.indBigDends(1:numel(axonsNew), 1) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, sprintf('dendrites_%s.mat', dendriteVersionNew));
    Util.saveStruct(saveFile, out);
    system(sprintf('chmod -w "%s"', saveFile));
end
