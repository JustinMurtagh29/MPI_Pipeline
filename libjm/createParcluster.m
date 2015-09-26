function cluster = createParcluster( clusterConfiguration, priority, memoryPerTask, jobStorageLocation)

    % Default values for all input arguments
    if nargin < 1
        clusterConfiguration = 'cpu';
    end
    if nargin < 2
        priority = -500;
    end
    if nargin < 3
        if strcmp(clusterConfiguration, 'cpu')
            memoryPerTask = 16;
        else
            memoryPerTask = 36;
        end
    end
    if nargin < 4
        jobStorageLocation = '/u/mberning/temp/';
    end

    % Generate generic cluster profile
    cluster = parallel.cluster.Generic();
    cluster.JobStorageLocation = jobStorageLocation;
    cluster.HasSharedFilesystem =  true;
    cluster.ClusterMatlabRoot = '/gaba/u/system/SLES11/soft/matlab/R2014bn/';
    cluster.OperatingSystem = 'unix';
    % Set all kind of utility functions
    cluster.GetJobStateFcn = @getJobStateFcn;
    cluster.DeleteJobFcn = @deleteJobFcn;
    % Set independent submit function parameter dependent on cluster configuration
    switch clusterConfiguration
        case 'm2090'
            cluster.IndependentSubmitFcn = {@independentSubmitFcn, ['-p ' num2str(priority) ' -pe openmp 1 -l h_vmem=' num2str(memoryPerTask) 'G,h_rt=10:00:00,num_m2090=1']};
        case 'k40'
            cluster.IndependentSubmitFcn = {@independentSubmitFcn, ['-p ' num2str(priority) ' -pe openmp 1 -l h_vmem=' num2str(memoryPerTask) 'G,h_rt=100:00:00,num_k40=1']};
        case 'gpu'
            cluster.IndependentSubmitFcn = {@independentSubmitFcn, ['-p ' num2str(priority) ' -pe openmp 1 -l h_vmem=' num2str(memoryPerTask) 'G,h_rt=100:00:00,num_k40=1']};
        case 'cpu'
            cluster.IndependentSubmitFcn = {@independentSubmitFcn, ['-p ' num2str(priority) ' -pe openmp 1 -l h_vmem=' num2str(memoryPerTask) 'G,h_rt=100:00:00']};
        otherwise
            error('getCluster:unknownClusterConfiguration', 'Unkown cluster configuration string passed');
    end
end

