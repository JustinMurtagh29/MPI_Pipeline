function job = calculateFeatures(p, voxelSource)
% Calculate SynEM features for SegEM classification output
%
% Written by
%   Manuel Berning <manuel.berning@brain.mpg.de>
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

% load fm
load([p.saveFolder 'SynapseClassifier.mat'], 'fm');
% needed changes
fm.areaT = 10;

% Construct data for job submission
fH = @connectEM.calculateFeaturesCube;
inputCell = arrayfun(@(x) {{x}}, 1:numel(p.local));

% Start feature calculation on gaba
cluster = Cluster.fetchCluster( ...
        'memory',24, ...
        'time','24:00:00');

job = Cluster.startJob( ...
    fH, inputCell, ...
    'name', strcat(lower(voxelSource), 'Features'), ...
    'sharedInputs', {p fm voxelSource}, ...
    'sharedInputsLocation', [1 3 4], ...
    'cluster', cluster);

end
