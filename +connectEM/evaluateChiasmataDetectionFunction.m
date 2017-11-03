function evaluateChiasmataDetectionFunction(aggloState,chiasmId,outputDir)
% This function loads chiasmata detection results and performs a simple
% quantitative analysis. These numbers should complement the impression
% gained from looking at random super-agglomerates (or parts thereof) with
% their detected chiasmata
% INPUT
% aggloState    string with the name of the aggloState that was used for
%               chiasmata detection
% chiasmId      string with id of the chiasmata detection run (datestring)
% outputDir     string with path to folder which will be created to safe
%               example nmls and results in
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de> / Marcel Beining

% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
% aggloState = 'dendrites_03_v2';
aggloFile = fullfile(rootDir, 'aggloState', [aggloState '.mat']);

% chiasmId = '20171019T155102';
% chiasmDir = '/tmpscratch/kboerg/chiasmata';
chiasmDir = fullfile(rootDir, '/aggloState/chiasmata/dendrites/round1/');
% chiasmataX20171017T044832_0
% outputDir = '/tmpscratch/mbeining/chiasma_03_v2';

%% load parameter (for skeleton)
param = fullfile(rootDir, 'allParameter.mat');
param = load(param);
param = param.p;

%% load all super-agglomerates
agglos = load(aggloFile);
if ~isempty(strfind(aggloFile,'axons'))
    agglos = agglos.axons(agglos.indBigAxons);
else
    agglos = agglos.dendrites(agglos.indBigDends);
end
aggloCount = numel(agglos);

%% load all chiasmata results
chiasmata = cell(aggloCount, 1);

tic;
for curIdx = 1:aggloCount
    curFile = fullfile( ...
        chiasmDir, ...
        sprintf('chiasmataX20171103T150231%s_%d', chiasmId, floor(curIdx / 100)), ...
        sprintf('visX20171103T150231%s_%d', chiasmId, curIdx), 'result.mat');
%         sprintf('chiasmataX%s_%d', chiasmId, floor(curIdx / 100)), ...
%         sprintf('visX%s_%d', chiasmId, curIdx), 'result.mat');
    
    curData = load(curFile);
    chiasmata{curIdx} = curData.output;
    
    if ~mod(curIdx, 500); fprintf('%d done\n', curIdx); end
end
toc;

%% start couting...
clearvars -except param outputDir chiasmata chiasmId aggloState;
mkdir(outputDir);
save(fullfile(outputDir,sprintf('chiasmata_%s_%s.mat',chiasmId,aggloState)))

% count number of nodes
nodeCount = cellfun(@(s) size(s.nodes, 1), chiasmata);
chiasmaNodeCounts = cellfun(@(s) sum(s.isIntersection), chiasmata);
chiasmaNodeFracs = chiasmaNodeCounts ./ nodeCount;

chiasmaCounts = cellfun(@(s) numel(s.ccCenterIdx), chiasmata);
chiasmaTotalCount = sum(chiasmaCounts);

chiasmaSizes = cell2mat(cellfun(@(s) ...
    cellfun(@numel, s.ccNodeIdx(:)), chiasmata, 'Uni', false));
chiasmaNrExits = cell2mat(cellfun(@(c) ...
    cellfun(@(p) size(p, 1), c.position), chiasmata, 'Uni', false));

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
temp.nrChiasmata = arrayfun(@(c) sum(chiasmaNrExits == c), temp.nrExits);
fprintf('\n'); disp(temp);

%% write out largest chiasmata and random chiasmata
[descVals, descIds] = sort(chiasmaSizes, 'descend');

% discard chiasmata with more than 500 nodes
% these cannot properly be inspected in webKNOSSOS anyway
descIds(descVals > 500) = [];
descVals(descVals > 500) = [];

% go back
invLUT = cat(1, 0, cumsum(chiasmaCounts));
randIdx = randperm(numel(descIds),10);
for curIdx = 1:10
    for n = 1:2
        switch n
            case 1
                curChiasmaIdx = descIds(curIdx);
            case 2
                curChiasmaIdx = descIds(randIdx(curIdx));
        end
        % convert back to indices
        curAggloIdx = find(invLUT >= curChiasmaIdx, 1) - 1;
        curChiasmaIdx = curChiasmaIdx - invLUT(curAggloIdx);
        
        curChiasma = chiasmata{curAggloIdx};
        curNodeIds = curChiasma.ccNodeIdx{curChiasmaIdx};
        if n == 1
            assert(numel(curNodeIds) == descVals(curIdx));
             curName = sprintf( ...
            'Superagglo_Largest_%d__Chiasma_%d', curAggloIdx, curChiasmaIdx);
        else
            assert(numel(curNodeIds) == descVals(randIdx(curIdx)));
            curName = sprintf( ...
            'Superagglo_Rand_%d__Chiasma_%d', curAggloIdx, curChiasmaIdx);
        end
        % build skeleton
        
        curFileName = fullfile(outputDir, strcat(lower(curName), '.nml'));
        
        skel = skeleton();
        skel = skel.addTree(curName, curChiasma.nodes, curChiasma.edges);
        skel = skel.addBranchpoint(curNodeIds);
        
        % show queries
        curPositions = curChiasma.position{curChiasmaIdx};
        curDirections = curChiasma.direction{curChiasmaIdx};
        
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
end
