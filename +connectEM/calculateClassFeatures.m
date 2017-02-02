function calculateClassFeatures( p )
% Calculate SynEM features for SegEM classification output

% load fm
load([p.saveFolder 'SynapseClassifier.mat'], 'fm');
% needed changes
fm.areaT = 10;
% not sure whether necessary
fm.voxelSize = [11.24 11.24 28];

% Construct data for job submission
fH = @connectEM.Seg.predictCube;
inputCell = cellfun(@(x){p, x, fm}, 'uni', 0);

% Start feature calculation on gaba
startCPU(fH, inputCell);

end
