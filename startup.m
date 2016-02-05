% Startup File pipeline repo
global CLUSTER_CPU CLUSTER_GPU;

% Add current directory to path
addpath(genpathGit(pwd));
% Set up gaba cluster configurations for usage
CLUSTER_CPU = Cluster.getCluster('-pe openmp 1', '-l h_vmem=32G', '-l h_rt=24:00:00');
CLUSTER_GPU = Cluster.getCluster('-pe openmp 1', '-l h_vmem=32G', '-l h_rt=24:00:00', 'num_k40=1');
