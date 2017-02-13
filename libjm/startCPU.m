function job = startCPU(fH, iC, jN, requiredMemory, group, rt);
    % Wrapper function for startJob.m used for backward compability
    if nargin < 6
	rt = 24;
    end
    if nargin < 4
        global CLUSTER_CPU;
    else
        CLUSTER_CPU = Cluster.getCluster('-pe openmp 1', ['-l h_vmem=' num2str(requiredMemory) 'G'], sprintf('-l s_rt=%02d:%02d:00 -l h_rt=%02d:%02d:30',floor(rt),round((rt-floor(rt))*60),floor(rt),round((rt-floor(rt))*60)), '-p -500');
    end
    if nargin < 5
        group = 1;
    end
    job = startJob(CLUSTER_CPU, fH, iC, jN, group);

end

