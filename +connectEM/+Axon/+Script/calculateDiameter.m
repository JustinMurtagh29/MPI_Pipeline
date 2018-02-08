% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

[outFile, outName] = fileparts(connFile);
outName = sprintf('%s_diameter.mat', outName);
outFile = fullfile(outFile, outName);
clear outName;

info = Util.runInfo();

%% loading data
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segSizes = Seg.Global.getSegToSizeMap(param);

segCentroids = Seg.Global.getSegToCentroidMap(param);
segCentroids = segCentroids .* param.raw.voxelSize;

segCov = load(fullfile(rootDir, 'globalSegmentPCA.mat'), 'covMat');
segCov = reshape(segCov.covMat, [], 3, 3);

%% calculate diameters
tic; fprintf('Calculating axon diameters... ');
[segIds, coms, diams] = Agglo.calculateDiameter( ...
    segSizes, segCentroids, segCov, conn.axons, 'nhoodThresh', 750);
fprintf('done!\n'); toc;

%% calculate node weights
% NOTE(amotta): To calculate the mean or median diameter of an axon, we
% need to know how much each measurement (i.e., node) contributes to the
% whole axon. Let's do this by calculating the sum of all half-edges for
% each node...
weights = cell(size(coms));
for curIdx = 1:numel(coms)
    curComs = coms{curIdx};
    
    if size(curComs, 1) == 1
        weights{curIdx} = 1;
    else
        curAdj = pdist(curComs);
        curAdj = sparse(squareform(curAdj));
        curAdj = graphminspantree(curAdj);
        
        % sanity check
        assert(size(curAdj, 1) == size(curComs, 1));
        
        curWeights = sum(curAdj, 1) + sum(curAdj, 2)';
        curWeights = full(curWeights(:)) / 2;
        
        weights{curIdx} = curWeights;
    end
end

%% save result
Util.save(outFile, info, segIds, coms, diams, weights);
