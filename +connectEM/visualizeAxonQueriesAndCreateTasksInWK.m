
    outputFolder = '/gaba/scratch/mberning/axonQueryGeneration/';
    % Load data again (id needed)
    %load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    %[graph, segmentMeta, borderMeta, globalSegmentPCA] = connectEM.loadAllSegmentationData(p);
    segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');

    % Load first batch for writing debug skeletons
    batch = load('/gaba/scratch/mberning/axonQueryGeneration/queriesMat/batch0001.mat');

    % Visualize axon queries as skeletons for debugging process
    idx = randperm(numel(batch.theseAxons), 100);
    connectEM.debugNewQueries(segmentMeta, batch.theseAxons, batch.q, [outputFolder 'queryVisualization/']);
    clear idx;

    % Write first batch to tasks on webKnossos
    % This needs the anaconda/2 module loaded before starting matlab to work (on gaba)
    for i=1:100
        batch = load(['/gaba/scratch/mberning/axonQueryGeneration/queriesMat/batch' num2str(i, '%.4i') '.mat']);
        tic;
        for j=1:length(batch.q.pos)
            command = sprintf('python ../wk-paper/RESCOPaaS/createTaskWithScript.py %u %u %u %.2f %.2f %.2f 1 MB_AxonFocusApiTest', ...
                batch.q.pos{j}(:), batch.q.angles{j}(:));
            system(command);
            pause(.1);
        end
        toc;
    end

