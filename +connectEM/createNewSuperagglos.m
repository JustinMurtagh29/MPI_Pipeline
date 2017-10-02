function createNewSuperagglos(param, state)
    % The key idea is that we want at most one flight path between any pair
    % of agglomerates during the generation of super-agglomerates.
    % 
    % In a first step we thus only look at flight paths which attach to
    % exactly two agglomerates. For each occuring pair of agglomerates we
    % want to pick exactly one flight path from one or multiple ones.
    % 
    % To do so, the flight paths connecting a given pair of agglomerates
    % are sorted by descending "reliability" and only the top one is
    % chosen. I've decided to use the number of involved endings as measure
    % for reliability:
    % 
    % Two is better than one is better than none.
    % 
    % What remains are the dangling flight paths: By definition these
    % attach to exactly one ending. To remove redundant flight paths we
    % allow at most one dangling flight path per ending.
    %
    % From an email from AM to the entire L4 team on Sunday, 01.10.2017.
    %
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    dataDir = fullfile(param.saveFolder, 'aggloState');

    [~, ~, suffix, axonVersion, axonVersionNew, casesToMerge] = ...
        connectEM.setQueryState(state);

    % Load current state of agglomerates
    agglos = load(fullfile( ...
        dataDir, sprintf('axons_%02d.mat', axonVersion)));
    superAgglos = agglos.axons(agglos.indBigAxons);
    
    % NOTE(amotta): `agglos` contains all super-agglomerates. `superAgglos`
    % is restricted to the large-enough ones.

    % Load linkages and cases for execusion
    load(fullfile(dataDir, sprintf('caseDistinctions%s.mat', suffix)), ...
        'linkagesAgglos', 'flightPaths', 'linkagesFlat');
    load(fullfile(dataDir, sprintf('attachedEndings%s.mat', suffix)), ...
        'flightsOfEndingCases', 'endingCaseDistinctions');
    
    % Choose cases for merging and generation of queries
    executedFlightPaths = flightsOfEndingCases( ...
        ismember(endingCaseDistinctions, casesToMerge));
    executedFlightPaths = unique(cat(2, executedFlightPaths{:})');
       
    % sanity checks
    assert(isequal(size(linkagesAgglos), size(linkagesFlat)));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.startAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.endAgglo));
    assert(size(linkagesAgglos, 1) == numel(flightPaths.ff.segIds));
    assert(~any(diff(structfun(@numel, flightPaths.ff))));
    assert(max(executedFlightPaths) < size(linkagesAgglos, 1));

    %% preparation for redundancy detection
    mask = false(size(linkagesAgglos, 1), 1);
    mask(executedFlightPaths) = true;

    % mask for flight that attach at both ends
    validAttach = ...
        @(ids) numel(ids) == 1 && ids(1) > 0;
    aggloMask = ...
        cellfun(validAttach, flightPaths.startAgglo) ...
      & cellfun(validAttach, flightPaths.endAgglo);
  
    % find ids of dangling and doubly attached flights
    dangIds = find(mask & ~aggloMask);
    aggloIds = find(mask & aggloMask);
    
    %% handle doubly attached flight paths
    aggloPairs = horzcat( ...
        cell2mat(flightPaths.startAgglo(aggloIds)), ...
        cell2mat(flightPaths.endAgglo(aggloIds)));
    aggloPairs = sort(aggloPairs, 2);

    % NOTE(amotta): Build a table with the pairs-of-agglomerates and pairs-
    % of-endings. By sorting this table DESCENDING ENDING IDS for each pair
    % of agglomerates and then taking the unique agglomerate pairs, we
    % favor ending attachments.
    endPairs = linkagesFlat(aggloIds, :);
    endPairs = sort(endPairs, 2);

    aggloEndPairs = horzcat(aggloPairs, endPairs);
   [aggloEndPairs, sortIds] = sortrows(aggloEndPairs, [1, 2, -3, -4]);
    aggloIds = aggloIds(sortIds);

   [uniAggloPairs, ~, pairFlightGroups] = ...
        unique(aggloEndPairs(:, 1:2), 'rows');
    pairFlightGroups = arrayfun( ...
        @(idx) aggloIds(pairFlightGroups == idx), ...
        (1:max(pairFlightGroups))', 'UniformOutput', false);

    %% handle dangling flight paths
    assert(all(linkagesFlat(dangIds, 1) > 0));
    assert(all(isnan(linkagesFlat(dangIds, 2))));
    dangEndIds = linkagesFlat(dangIds, 1);

   [~, ~, dangFlightGroups] = unique(dangEndIds);
    dangFlightGroups = arrayfun( ...
        @(idx) dangIds(dangFlightGroups == idx), ...
        (1:max(dangFlightGroups))', 'UniformOutput', false);

    %% put it all together
    flightGroups = cat(1, pairFlightGroups, dangFlightGroups);
    uniFlightIds = cellfun(@(ids) ids(1), flightGroups);

   [uniFlightIds, sortIds] = sort(uniFlightIds);
    flightGroups = flightGroups(sortIds);

    eqClasses = Graph.findConnectedComponents(uniAggloPairs, true, true);
    eqClasses = [eqClasses; num2cell(setdiff( ...
        1:numel(superAgglos), cell2mat(eqClasses)))'];

    %% build super-agglomerates
    flightPaths.startAgglo = ...
        flightPaths.startAgglo(uniFlightIds);
    flightPaths.endAgglo = ...
        flightPaths.endAgglo(uniFlightIds);
    flightPaths.ff = structfun( ...
        @(x) x(uniFlightIds), flightPaths.ff, ...
        'UniformOutput', false);

    axonsNew = connectEM.mergeSuperagglosBasedOnFlightPath( ...
        superAgglos, eqClasses, flightPaths.startAgglo, ...
        flightPaths.endAgglo, flightPaths.ff);

    % sanity check
    assert(all(arrayfun(@(x) ...
        numel(Graph.findConnectedComponents(x.edges)) == 1, ...
        axonsNew(arrayfun(@(x) size(x.nodes, 1), axonsNew) > 1))));

    out = struct;
    out.param = param;
    out.state = state;
    out.gitInfo = Util.gitInfo();
    
    out.eqClasses = eqClasses;
    out.flightGroups = flightGroups;
    
    out.axons = cat(1, axonsNew, agglos.axons(~agglos.indBigAxons));
    out.indBigAxons = false(length(out.axons),1);
    out.indBigAxons(1:numel(axonsNew), 1) = true;

    % Save super agglos and deprive writing permission
    saveFile = fullfile(dataDir, sprintf('axons_%s.mat', axonVersionNew));
    Util.saveStruct(saveFile, out);
    system(sprintf('chmod -w "%s"', saveFile));
end
