function job = detectChiasmataSuperSuper( ...
    p, inputFile, outputDir, useSphereClustering)

if ~exist('useSphereClustering', 'var')
    % set to true for alternate approach of clustering on sphere
    % (detectChiasmata vs. detectChiasmataSphereClustering)
    useSphereClustering = false;
end

%for a clean version of findCCaccordingToGraph
addpath('/gaba/u/kboerg/code/manuelCode/games');
functionH = @connectEM.detectChiasmataSuper;
inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);

% set id for detected chiasmata
p.inputFile = inputFile;
p.outputDir = outputDir;
p.chiasmataVersion = datestr(now, 30);

% Some parameter for algorithm
p.sphereRadiusOuter = 10000; % in nm
p.sphereRadiusInner = 1000; % in nm
p.minNodeDist = 2000; % in nm
p.clusterSize = 2000; % in nm
p.voxelSize = p.raw.voxelSize;

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata', ...
    'sharedInputs', {p, useSphereClustering}, ...
    'sharedInputsLocation', [2, 3], ...
    'cluster', cluster);
end
