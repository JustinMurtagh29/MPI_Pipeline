
    % Load data again (id needed)
    %load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    %[graph, segmentMeta, borderMeta, globalSegmentPCA] = connectEM.loadAllSegmentationData(p);

    % Load first batch for writing debug skeletons
    batch = load('/gaba/scratch/mberning/axonQueryGeneration/queriesMat/batch0001.mat');

    % Visualize axon queries as skeletons for debugging process
    idx = randperm(numel(theseAxons), 100);
    connectEM.debugNewQueries(segmentMeta, batch.theseAxons, batch.q, [outputFolder 'queryVisualization/']);
    clear idx;

    % Start working on using RESCOPaaS like access to webKnossos REST
    % K I guess we need to start working on wklib again to transfer q.pos + q.angles into new tasks (with scripts)
