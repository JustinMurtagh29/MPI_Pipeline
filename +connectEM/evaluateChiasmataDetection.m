% This script loads chiasmata detection results and performs a simple
% quantitative analysis. These numbers should complement the impression
% gained from looking at random super-agglomerates (or parts thereof) with
% their detected chiasmata
% 
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_10_a.mat');

chiasmId = '20171104T131021';
chiasmDir = '/tmpscratch/amotta/l4/chiasma-detection';

outputDir = '/home/amotta/Desktop';

%% load parameter (for skeleton)
param = fullfile(rootDir, 'allParameter.mat');
param = load(param);
param = param.p;

%% load all super-agglomerates
agglos = load(axonFile);
agglos = agglos.axons(agglos.indBigAxons);
aggloCount = numel(agglos);

%% load all chiasmata results
chiasmata = ...
    connectEM.Chiasma.Detect.collectResults( ...
        chiasmDir, chiasmId, aggloCount);
connectEM.Chiasma.Detect.evaluate(chiasmata, agglos);

%% write out largest chiasmata
mkdir(outputDir);
[descVals, descIds] = sort(chiasmaSizes, 'descend');

% discard chiasmata with more than 500 nodes
% these cannot properly be inspected in webKNOSSOS anyway
descIds(descVals > 500) = [];
descVals(descVals > 500) = [];

% go back
invLUT = cat(1, 0, cumsum(chiasmaCounts));

for curIdx = 1:10
    curChiasmaIdx = descIds(curIdx);
    
    % convert back to indives
    curAggloIdx = find(invLUT >= curChiasmaIdx, 1) - 1;
    curChiasmaIdx = curChiasmaIdx - invLUT(curAggloIdx);
    
    curAgglo = chiasmata{curAggloIdx};
    curNodeIds = curAgglo.ccNodeIdx{curChiasmaIdx};
    assert(numel(curNodeIds) == descVals(curIdx));
    
    % build skeleton
    curName = sprintf( ...
        'Superagglo_%d__Chiasma_%d', curAggloIdx, curChiasmaIdx);
    curFileName = fullfile(outputDir, strcat(lower(curName), '.nml'));
    
    skel = skeleton();
    skel = skel.addTree(curName, curAgglo.nodes, curAgglo.edges);
    skel = skel.addBranchpoint(curNodeIds);
    
    % show queries
    curPositions = curAgglo.position{curChiasmaIdx};
    curDirections = curAgglo.direction{curChiasmaIdx};
    
    for curQueryIdx = 1:size(curPositions, 1)
        curStartPos = curPositions(curQueryIdx, :);
        curEndPos = curStartPos - curDirections(curQueryIdx, :);
        
        skel = skel.addTree( ...
            sprintf( ...
                'Chiasma_%d__Query_%d', ...
                curChiasmaIdx, curQueryIdx), ...
            cat(1, curStartPos, curEndPos));
    end
    
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel.write(curFileName);
end
