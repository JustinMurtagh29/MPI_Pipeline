function [job1, job2] = testClusterReproduceability()

    functionH = @someRandomMatrixMultiplicationAndSumming;
    % Start once on an "oldgaba" node
    clusterCPU = Cluster.getCluster( ...
        '-pe openmp 1', ...
        ['-p ' num2str(0)], ...
        '-l hostname=gaba001');
    job1 = Cluster.startJob(functionH, {{}}, 'cluster', clusterCPU, 'name', 'reproduce', 'numOutputs', 3);
    % Start once on an "newgaba" node
    clusterCPU = Cluster.getCluster( ...
        '-pe openmp 1', ...
        ['-p ' num2str(0)], ...
        '-l hostname=gaba024');
    job2 = Cluster.startJob(functionH, {{}}, 'cluster', clusterCPU, 'name', 'reproduce', 'numOutputs', 3);
    % Wait until both jobs have finished
    wait(job1, 'finished');
    wait(job2, 'finished');
    % Collect outputs
    result1 = getAllOutputArguments(job1);
    result2 = getAllOutputArguments(job2);
    % Compare outputs
    display(all(result1{1}(:) == result2{1}(:)));
    display(all(result1{2}(:) == result2{2}(:)));
    display(all(result1{3}(:) == result2{3}(:)));
end

function [matrix, matrix10, matrixSum] = someRandomMatrixMultiplicationAndSumming()
    rng default;
    matrix = rand(100, 1000, 1000, 'single').*1e8;
    matrix10 = matrix.^10;
    matrixSum = sum(matrix10(:));
end

