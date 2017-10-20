function [job,chiasmataVersion] = detectChiasmataSuperSuper(p, inputFile, useSphereClustering,visualization)

if ~exist('useSphereClustering', 'var')
    % set to true for alternate approach of clustering on sphere
    % (detectChiasmata vs. detectChiasmataSphereClustering)
    useSphereClustering = false;
end

%for a clean version of findCCaccordingToGraph
addpath('/gaba/u/kboerg/code/manuelCode/games');
functionH = @connectEM.detectChiasmataSuper;
inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);

% Some parameter for algorithm
if ~isfield(p,'sphereRadiusOuter')
    p.sphereRadiusOuter = 5000; % in nm
end
if ~isfield(p,'sphereRadiusInner')
    p.sphereRadiusInner = 2000; % in nm
end
if ~isfield(p,'minDistNode')
    p.minDistNode = 3000; % in nm
end
fprintf('OuterSphere: %d nm ; InnerSphere: %d nm ; minDistNode: %d nm\n',p.sphereRadiusOuter, p.sphereRadiusInner, p.minDistNode);
p.voxelSize = p.raw.voxelSize;
% set id for detected chiasmata
p.inputFile = inputFile;
chiasmataVersion = datestr(now, 30);
p.chiasmataVersion = chiasmataVersion;
cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata', ...
    'sharedInputs', {p, useSphereClustering,visualization}, ...
    'sharedInputsLocation', [2, 3, 4], ...
    'cluster', cluster);
