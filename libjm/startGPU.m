function job = startGPU(fH, iC, jN);
    % Wrapper function for startJob.m used for backward compability

    global CLUSTER_GPU;
    job = startJob(CLUSTER_GPU, fH, iC, jN);

end

