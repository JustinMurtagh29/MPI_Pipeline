function flightEndingOverlapRun(param)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

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
    
    % load flight paths
    m = load(fullfile(dataDir, 'AxonFlightPaths.mat'), 'ff');
    flightNodes = m.ff.nodes;
    
    % load results from flight path evaluation
    m = load(fullfile(dataDir, 'AxonQueryOverlaps.mat'), 'results');
    flightResults = m.results;
    clear m;

    %% run main function
    doIt = @(ids) connectEM.flightEndingOverlap( ...
            param, origAgglos, endings, flightNodes, ids, superAgglos);
        
    out = struct;
    out.startEndingOverlaps = doIt(flightResults.startAgglo);
    out.endEndingOverlaps = doIt(flightResults.endAgglos);
    
    %% save result
    outFile = fullfile(dataDir, 'AxonEndingOverlaps.mat');
    Util.saveStruct(outFile, out);
end
