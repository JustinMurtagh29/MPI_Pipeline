function job = calculateSegFeatures( p, fmloc )
% Calculate SynEM features for SegEM output

% load fm
if nargin < 2
    load(fullfile(p.saveFolder,'SegmentFeatureMap.mat'),'fm')
else
    load(fmloc, 'fm');
emd
% needed changes
% fm.areaT = 200;  % segment voxel size should exceed 200 voxels
% fm.numSubvolumes = 0; % only each segment itself

% Construct data for job submission
fH = @alignEM.calculateSegFeaturesCube;
inputCell = cellfun(@(x){x}, num2cell(1:numel(p.local),1), 'uni', 0);

% Start feature calculation on gaba
cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p -500', ...
        '-l h_vmem=24G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
job = Cluster.startJob(fH, inputCell, 'name', 'segFeatures', 'sharedInputs', {p fm}, 'sharedInputsLocation', [1 3], 'cluster', cluster);

end
