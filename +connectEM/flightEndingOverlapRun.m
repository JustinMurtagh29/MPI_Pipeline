function flightEndingOverlapRun(param,state)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% Set current state of queries
    [skeletonFolders, suffixFlightPaths, suffix] = connectEM.setQueryState(state);    

    %% load all the input data

    dataDir = fullfile(param.saveFolder, 'aggloState');

    % load axon super-agglomerates
    m = load(fullfile(dataDir, 'axons_04.mat'));

    % IMPORTANT(amotta): These must be exactly the agglomerates to which
    % `startAgglo` and `endAgglo` of `flightResults` refers to!
    superAgglos = m.axons(m.indBigAxons);

    % IMPORTANT(amotta): These must be exactly the agglomerates to which
    % `endings.aggloIds` refers to!
    origAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, m.axons, 'UniformOutput', false);

    % load endings
    endings = load(fullfile(dataDir, 'axonEndings.mat'));

    m = load(fullfile(dataDir, strcat('axonPostQueryAnalysisState',suffix,'.mat')), 'ff', 'results');
    % load flight paths
    flightNodes = m.ff.nodes;
    % load results from flight path evaluation
    flightResults = m.results;
    clear m;

    %% run main function
    doIt = @(ids) connectEM.flightEndingOverlap( ...
            param, origAgglos, endings, flightNodes, ids, superAgglos);

    out = struct;
    out.startEndingOverlaps = doIt(flightResults.startAgglo);
    out.endEndingOverlaps = doIt(flightResults.endAgglo);

    %% save result
    outFile = fullfile(dataDir, strcat('axonEndingOverlaps',suffix,'.mat'));
    Util.saveStruct(outFile, out);

end
