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
    m = load(fullfile(dataDir, 'axonFlightPaths.mat'), 'ff');
    m.ff = structfun(@(x)x(cellfun(@isempty, m.ff.comments)), m.ff, 'uni', 0);
    m.ff = structfun(@(x)x(~cellfun(@isempty, m.ff.startNode)), m.ff, 'uni', 0);
    flightNodes = m.ff.nodes;

    % load results from flight path evaluation
    m = load(fullfile(dataDir, 'axonQueryOverlaps.mat'), 'results');
    flightResults = m.results;
    clear m;

    %% run main function
    doIt = @(ids) connectEM.flightEndingOverlap( ...
            param, origAgglos, endings, flightNodes, ids, superAgglos);

    out = struct;
    out.startEndingOverlaps = doIt(flightResults.startAgglo);
    out.endEndingOverlaps = doIt(flightResults.endAgglo);

    %% save result
    outFile = fullfile(dataDir, 'axonEndingOverlaps.mat');
    Util.saveStruct(outFile, out);


    %% Display some statistics
    % Endings statistics:
    endingClusters = endings.borderClusters;
    clusterSizes = cellfun(@max, endingClusters);
    singleEnding = sum(clusterSizes == 1);
    display([num2str(singleEnding./numel(clusterSizes)*100, '%.2f') '% of agglomerates have just one single ending']);
    display([num2str(singleEnding) ' in total']);

    startEndings = unique(cell2mat(out.startEndingOverlaps));
    endEndings = unique(cell2mat(out.endEndingOverlaps));
    totalEndings = union(startEndings,endEndings);

    display([num2str(numel(startEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached at start']);
    display([num2str(numel(startEndings)) ' in total']);
    display([num2str(numel(endEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached at end']);
    display([num2str(numel(endEndings)) ' in total']);

    display([num2str(numel(totalEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached']);
    display([num2str(numel(totalEndings)) ' in total']);

    % Flight path statistics:
    display([num2str(sum(~cellfun('isempty',out.startEndingOverlaps))./numel(flightEndingOverlap.starts)*100, '%.2f')...
        '% of flight paths attach at start']);
    display([num2str(sum(~cellfun('isempty',out.startEndingOverlaps))) ' in total']);
    display([num2str(sum(~cellfun('isempty',out.endEndingOverlaps))./numel(flightEndingOverlap.starts)*100, '%.2f')...
        '% of flight paths attach at end']);
    display([num2str(sum(~cellfun('isempty',out.endEndingOverlaps))) ' in total']);

    totalEndingOverlaps = cellfun( ...
        @union, out.startEndingOverlaps, ...
        out.endEndingOverlaps, 'UniformOutput', false);

    display([num2str(sum(~cellfun('isempty',totalEndingOverlaps))./numel(totalEndingOverlaps)*100, '%.2f') '% of flight paths attach at ending']);
    display([num2str(sum(~cellfun('isempty',totalEndingOverlaps))) ' in total']);

end
