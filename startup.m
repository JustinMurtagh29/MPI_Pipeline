% Add current directory and its sub-directories to path
addpath(genpathGit(fileparts(mfilename('fullpath'))));

% Set up gaba cluster configurations for usage
% 12 GB normal job limit, 24 hour run time and -500 priority
global CLUSTER_CPU;
CLUSTER_CPU = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p -500', ...
    '-l h_vmem=12G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');

% If asking for a GPU, set normal (0) priority
global CLUSTER_GPU;
CLUSTER_GPU = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=36G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00', ...
    '-l num_k40=1');

% Mark as ready
global PIPELINE_READY;
PIPELINE_READY = true;


% By default dataset is not a retina so p.retina = false. REmoce retina directories from the matlab path
global p_retina;
p_retina = false;

if p.retina
rmpath([pwd '/CNN']);
rmpath([pwd '/segmentation']);
else
rmpath([pwd '/retinaCNN']);
rmpath([pwd '/retinaSegmentation']);
end

