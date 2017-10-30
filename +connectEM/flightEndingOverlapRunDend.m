function flightEndingOverlapRunDend(param,state)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% Set current state of queries
    [~, suffix, dendriteVersion] = connectEM.setQueryState(state);    

    %% load all the input data

    dataDir = fullfile(param.saveFolder, 'aggloState');

    % load axon super-agglomerates
    m = load(fullfile(dataDir, strcat('dendrites_',dendriteVersion)));

    % IMPORTANT(amotta): These must be exactly the agglomerates to which
    % `startAgglo` and `endAgglo` of `flightResults` refers to!
    superAgglos = m.dendrites(m.indBigDends);

    % IMPORTANT(amotta): These must be exactly the agglomerates to which
    % `endings.aggloIds` refers to!
    origAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, m.dendrites, 'UniformOutput', false);

    % load endings
    endings = load(fullfile(dataDir, strcat('dendriteEndings_',dendriteVersion,'.mat')));

    m = load(fullfile(dataDir, strcat('dendritePostQueryAnalysisState_',suffix,'.mat')), 'ff', 'results');
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

    %% save result and deprive writing permission
    outFile = fullfile(dataDir, strcat('dendriteEndingOverlaps_',suffix,'.mat'));
    Util.saveStruct(outFile, out);
    system(['chmod -w ' outFile])

end
