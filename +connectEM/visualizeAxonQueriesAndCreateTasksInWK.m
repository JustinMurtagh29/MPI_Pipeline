    % Visualize axon queries as skeletons for debugging process
    idx = randperm(numel(axonsNew), 100);
    connectEM.debugNewQueries(segmentMeta, axonsNew(idx), q(idx), [outputFolder 'queryVisualization/']);
    clear idx;

    % Start working on using RESCOPaaS like access to webKnossos REST

