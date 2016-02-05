% Startup File pipeline repo
global CLUSTER_CPU CLUSTER_GPU;

% Add current directory to path
addpath(genpathGit(pwd));
% Set up gaba cluster configurations for usage
% 12 GB normal job limit, 24 hour run time and -500 priority
CLUSTER_CPU = Cluster.getCluster('-pe openmp 1', '-l h_vmem=12G', '-l h_rt=24:00:00', '-p -500');
% If asking for a GPU, set normal (0) priority
CLUSTER_GPU = Cluster.getCluster('-pe openmp 1', '-l h_vmem=36G', '-l h_rt=24:00:00', '-l num_k40=1');

