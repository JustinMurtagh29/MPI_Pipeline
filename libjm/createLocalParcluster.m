function cluster = createLocalParcluster( jobStorageLocation, numWorkers )

    % Default values for all input arguments
    if nargin < 1
        jobStorageLocation = '/gaba/scratch/mberning/matlab-jobs/';
    end
    if nargin < 2
       numWorkers = 4;
    end

    % Create directory to store job data if it does not exist
    if ~exist(jobStorageLocation, 'dir')
        mkdir(jobStorageLocation);
    end
    
    % Generate local cluster profile
    cluster = parallel.cluster.Local();
    cluster.NumWorkers = numWorkers;
    cluster.JobStorageLocation = jobStorageLocation;

end

