function job = startCPU(fH, iC, jN);
    % Wrapper function for startJob.m used for backward compability

    global CLUSTER_CPU;
    job = startJob(CLUSTER_CPU, fH, iC, jN);

end

