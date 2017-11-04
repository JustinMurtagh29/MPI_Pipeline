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
chiasmata = cell(aggloCount, 1);

tic;
for curIdx = 1:aggloCount
    curFile = fullfile( ...
        chiasmDir, ...
        sprintf('chiasmataX%s_%d', chiasmId, floor(curIdx / 100)), ...
        sprintf('visX%s_%d', chiasmId, curIdx), 'result.mat');
    
    curData = load(curFile);
    chiasmata{curIdx} = curData.output;
    
    if ~mod(curIdx, 500); fprintf('%d done\n', curIdx); end
end
toc;

%% start couting...
clearvars -except param agglos outputDir chiasmata;

% count number of nodes
nodeCount = cellfun(@(s) size(s.nodes, 1), chiasmata);
chiasmaNodeCounts = cellfun(@(s) sum(s.isIntersection), chiasmata);
chiasmaNodeFracs = chiasmaNodeCounts ./ nodeCount;

chiasmaCounts = cellfun(@(s) numel(s.ccCenterIdx), chiasmata);
chiasmaTotalCount = sum(chiasmaCounts);

% TODO(amotta): Use table instead
chiasmaSizes = cell2mat(cellfun(@(s) ...
    cellfun(@numel, s.ccNodeIdx(:)), chiasmata, 'Uni', false));
chiasmaNrExits = cell2mat(cellfun(@(c) ...
    cellfun(@(p) size(p, 1), c.position), chiasmata, 'Uni', false));

if isfield(agglos, 'solvedChiasma')
    chiasmaSolved = cell2mat(arrayfun(@(c, a) arrayfun( ...
        @(idx) double(a.solvedChiasma(idx)), c{1}.ccCenterIdx(:)), ...
        chiasmata, agglos, 'Uni', false));
    chiasmaSolved = logical(chiasmaSolved);
else
    chiasmaSolved = false(size(chiasmaNrExits));
end

chiasmaAggloMask = chiasmaCounts > 0;
chiasmaAggloFrac = mean(chiasmaAggloMask);
chiasmaAggloTotal = sum(chiasmaAggloMask);

%% report
fprintf('\nTotal number of chiasmata: %d\n', chiasmaTotalCount);

temp = sort(chiasmaCounts, 'descend');
fprintf('\nLargest number of chiasmata per agglomerate:\n');
fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

temp = sort(chiasmaSizes, 'descend');
fprintf('\nLargest chiasmata (in number of nodes):\n');
fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

temp = sort(chiasmaNodeCounts, 'descend');
fprintf('\nLargest number of chiasma nodes per agglomerate:\n');
fprintf('%s\n', strjoin(arrayfun(@num2str, temp(1:10), 'Uni', false), ', '));

temp = table;
temp.nrExits = (4:max(chiasmaNrExits))';
temp.nrChiasmata = arrayfun(@(c) ...
    sum(chiasmaNrExits == c), temp.nrExits);
temp.nrMarkedSolved = arrayfun(@(c) ...
    sum(chiasmaSolved(chiasmaNrExits == c)), temp.nrExits);
fprintf('\n'); disp(temp);

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
