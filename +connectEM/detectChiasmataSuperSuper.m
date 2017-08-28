function detectChiasmataSuperSuper(p)

addpath('/gaba/u/kboerg/code/manuelCode/games') %for a clean version of findCCaccordingToGraph
functionH = @connectEM.detectChiasmataSuper;
inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata', ...
    'sharedInputs', {p},  'sharedInputsLocation', 2, ...
    'cluster', cluster);
