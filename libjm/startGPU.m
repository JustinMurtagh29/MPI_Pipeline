function startGPU(fH, iC, jN);
    % Wrapper function for startJob.m used for backward compability

    global clusterGPU;
    startJob(clusterGPU, fH, iC, jN);
end

