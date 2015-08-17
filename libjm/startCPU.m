function startCPU(fH, iC, jN);
    % Wrapper function for startJob.m used for backward compability

    global clusterCPU;
    startJob(clusterCPU, fH, iC, jN);
end

