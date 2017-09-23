function createNewSuperagglos(param,state)
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    
    dataDir = fullfile(param.saveFolder, 'aggloState');

    [~, ~, suffix, axonVersion, axonVersionNew, casesToMerge] = connectEM.setQueryState(state);

    % Load current state of agglomerates
    agglos = load(fullfile(dataDir, strcat('axons_',num2str(axonVersion,'%.2i'),'.mat')));
    superAgglos = agglos.axons(agglos.indBigAxons);
    
    % NOTE(amotta): `agglos` contains all super-agglomerates. `superAgglos`
    % is restricted to the large-enough ones.

    % Load linkages and cases for execusion
    load(fullfile(dataDir, strcat('caseDistinctions',suffix,'.mat')),...
        'linkagesAgglos','flightPaths','linkagesFlat');
    load(fullfile(dataDir, strcat('attachedEndings',suffix,'.mat')),'flightsOfEndingCases','endingCaseDistinctions');
    
    % Choose cases for merging and generation of queries
    if isempty(casesToMerge)
        casesToMerge = [1:6, 8:14];
    end
    executedFlightPaths = flightsOfEndingCases(ismember(endingCaseDistinctions, casesToMerge));
    executedFlightPaths = unique(cat(2,executedFlightPaths{:})');
       
    % sanity checks
    assert(isequal(size(linkagesAgglos), size(linkagesFlat)));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.endAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
    assert(~any(diff(structfun(@numel, flightPaths.ff))));
    assert(max(executedFlightPaths) < size(linkagesAgglos, 1));
    
    execusion = false(size(linkagesAgglos,1),1);
    execusion(executedFlightPaths) = true;
    
    executedFlightPaths = execusion;
    clear execusion;

    % NOTE(amotta):
    % Q: `superAgglos` contains only the large-enough axons. Is it also
    %    true that the agglomerate IDs in `linkagesAgglos` refer to the
    %    large-enough axons?
    % A: The `linkagesAgglos` are directly generated from `startAgglo` and
    %    `endAgglo` of `flightPaths`. It seems like these fields also only
    %    refer to the large-enough axons.
    % See:
    %    https://gitlab.mpcdf.mpg.de/connectomics/pipeline/blob/master/+connectEM/getAggloQueryOverlapB.m#L18

    % Find connected components of equivalent classes
    % NOTE(amotta): Only use flights which attach at both end to
    % agglomerates. This makes sense.
    edgesCC = linkagesAgglos(executedFlightPaths,:);
    edgesCC(any(isnan(edgesCC), 2), :) = [];
    edgesCC = sort(edgesCC, 2);
    
    eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
    eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(superAgglos), cell2mat(eqClassCC)))'];

    % NOTE(amotta): Add flight path ID for tracking
    flightPaths.id = reshape(1:numel(flightPaths.startAgglo), [], 1);
    
    % Use only those flight paths for execusion with the chosen case IDs
    % NOTE(amotta): Dimensions agree!
    flightPaths.id(~executedFlightPaths) = [];
    flightPaths.startAgglo(~executedFlightPaths,:) = [];
    flightPaths.endAgglo(~executedFlightPaths,:) = [];
    flightPaths.ff = structfun(@(x)x(executedFlightPaths), flightPaths.ff, 'uni', 0);
    linkagesFlat(~executedFlightPaths, :) = [];

    % Eliminate duplicate flight paths
    % NOTE(amotta): Flight paths do **not** require to attach to
    % agglomerates at both ends. This is true for dangling flight queries,
    % for example.
    [~, positionUniques] = unique(sort(linkagesFlat,2), 'rows');
    duplicates = true(size(linkagesFlat,1),1);
    duplicates(positionUniques) = 0;
    
    % NOTE(amotta): No-endings (i.e., ending #0) cannot be compared for
    % equality. Ending zero of one flight path **is not the same** as
    % the ending zero of another flight path. These are the flight paths
    % producing case E5c endings.
    duplicates(any(linkagesFlat == 0, 2)) = false;
    
    % NOTE(amotta): Dimensions agree!
    flightPaths.id(duplicates) = [];
    flightPaths.startAgglo(duplicates,:) = [];
    flightPaths.endAgglo(duplicates,:) = [];
    flightPaths.ff = structfun(@(x)x(~duplicates), flightPaths.ff, 'uni', 0);
    
    % NOTE(amotta): At this point each agglomerate should still be a
    % single connected component. That's what Marcel checked. Can confirm.
    
    % sanity check
    % NOTE(amotta): passes!
    % assert(all(arrayfun(@(x) ...
    %     numel(Graph.findConnectedComponents(x.edges)) == 1, ...
    %     superAgglos(arrayfun(@(x) size(x.nodes, 1), superAgglos) > 1))));

    % Merging to super agglos
    axonsNew = connectEM.mergeSuperagglosBasedOnFlightPath( ...
        superAgglos, eqClassCCfull, flightPaths.startAgglo, ...
        flightPaths.endAgglo, flightPaths.ff);
    
    % sanity check
    assert(all(arrayfun(@(x) ...
        numel(Graph.findConnectedComponents(x.edges)) == 1, ...
        axonsNew(arrayfun(@(x) size(x.nodes, 1), axonsNew) > 1))));
    
    % Concatenate small axons below 5 um
    axons = cat(1,axonsNew, agglos.axons(~agglos.indBigAxons));
    indBigAxons = false(length(axons),1);
    indBigAxons(1:length(axonsNew),1) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, strcat('axons_',axonVersionNew,'.mat'));
    save(saveFile, 'axons','indBigAxons');
    system(['chmod -w ' saveFile]);
end
