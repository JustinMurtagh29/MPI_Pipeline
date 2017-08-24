function flightEndingOverlapRun(param)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% load all the input data
    dataDir = fullfile(param.saveFolder, 'aggloState');

    % load axon super-agglomerates
    m = load(fullfile(dataDir, 'axons_04.mat'));
    superAgglos = m.axons;

    % NOTE(amotta): In this particular case, the "original agglomerates"
    % used for the ending detection are identical to the agglomerates on
    % which flight paths might attach.
    origAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, superAgglos, 'UniformOutput', false);

    % load endings
    endings = load(fullfile(dataDir, 'axonEndings.mat'));

    m = load(fullfile(dataDir, 'AxonFlightPaths.mat'), 'ff');
    flightNodes = m.ff.nodes;

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
