function flightEndingOverlap = flightEndingOverlapRun()

    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% load all the input data
    m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    param = m.p;

    m = load(fullfile(param.saveFolder, 'aggloState/', 'axons_04.mat'));
    origAgglos = m.axons;

    endings = load(fullfile(param.saveFolder, 'aggloState/', 'axonEndings.mat'));

    m = load(fullfile(param.saveFolder, 'aggloState/', 'AxonFlightPaths.mat'), 'ff');
    flights = m.ff;

    m = load(fullfile(param.saveFolder, 'aggloState/', 'AxonQueryOverlaps.mat'), 'results')
    flightResults = m.results;

    m = load(fullfile(param.saveFolder, 'aggloState/', 'axons_04.mat'));
    superAgglos = m.axons;

    clear m;

    %% run main function
    flightNodes = flights.nodes;
    
    flightEndingOverlap.starts = connectEM.flightEndingOverlap( ...
            param, origAgglos, endings, flightNodes, flightResults.startAgglo, superAgglos);

    flightEndingOverlap.ends = connectEM.flightEndingOverlap( ...
            param, origAgglos, endings, flightNodes, flightResults.endAgglo, superAgglos);
        
    %% Display some statistics
    % Endings statistics:
    endingClusters = endings.borderClusters;
    clusterSizes = cellfun(@max, endingClusters);
    singleEnding = sum(clusterSizes == 1);
    display([num2str(singleEnding./numel(clusterSizes)*100, '%.2f') '% of agglomerates have just one single ending']);
    display([num2str(singleEnding) ' in total']);
    
    startEndings = unique(cell2mat(flightEndingOverlap.starts));
    endEndings = unique(cell2mat(flightEndingOverlap.ends));
    totalEndings = union(startEndings,endEndings);
    
    display([num2str(numel(startEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached at start']);
    display([num2str(numel(startEndings)) ' in total']);
    display([num2str(numel(endEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached at end']);
    display([num2str(numel(endEndings)) ' in total']);
    
    display([num2str(numel(totalEndings)./numel(endingClusters)*100, '%.2f') '% of endings have flight path attached']);
    display([num2str(numel(totalEndings)) ' in total']);
    
    % Flight path statistics:
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.starts))./numel(flightEndingOverlap.starts)*100, '%.2f')...
        '% of flight paths attach at start']);
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.starts))) ' in total']);
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.ends))./numel(flightEndingOverlap.starts)*100, '%.2f')...
        '% of flight paths attach at end']);
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.ends))) ' in total']);
    
    flightEndingOverlap.total = cellfun( ...
        @union, flightEndingOverlap.starts, ...
        flightEndingOverlap.ends, 'UniformOutput', false);
    
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.total))./numel(flightEndingOverlap.total)*100, '%.2f') '% of flight paths attach at ending']);
    display([num2str(sum(~cellfun('isempty',flightEndingOverlap.total))) ' in total']);
    
    
    

