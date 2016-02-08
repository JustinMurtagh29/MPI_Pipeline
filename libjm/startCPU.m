function job = startCPU(fH, iC, jN, requiredMemory);
    % Wrapper function for startJob.m used for backward compability

    if nargin < 4
        global CLUSTER_CPU;
    else
        CLUSTER_CPU = Cluster.getCluster('-pe openmp 1', ['-l h_vmem=' num2str(requiredMemory) 'G'], '-l h_rt=24:00:00', '-p -500');
    end
    job = startJob(CLUSTER_CPU, fH, iC, jN);

end

