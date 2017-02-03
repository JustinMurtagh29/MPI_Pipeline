function calculateClassFeatures( p )
% Calculate SynEM features for SegEM classification output

% load fm
load([p.saveFolder 'SynapseClassifier.mat'], 'fm');
% needed changes
fm.areaT = 10;

% Construct data for job submission
fH = @connectEM.calculateClassFeaturesCube;
inputCell = cellfun(@(x){x}, mat2cell(1:numel(p.local),1,ones(1,numel(p.local))), 'uni', 0);

% Start feature calculation on gaba
cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p -500', ...
        '-l h_vmem=24G', ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');
Cluster.startJob(fH, inputCell, 'name', 'classFeatures', 'sharedInputs', {p fm}, 'sharedInputsLocation', [1 3], 'cluster', cluster);

end
