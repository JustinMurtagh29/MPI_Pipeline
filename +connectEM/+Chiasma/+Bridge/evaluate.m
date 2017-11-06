% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
chiasmataFile = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a', ...
    '20171104T184018_chiasmata.mat');

outputDir = '/home/amotta/Desktop/random-chiasmata';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

data = load(chiasmataFile);
axonFile = data.info.param.axonFile;
chiasmaParam = data.info.param.chiasmaParam;
chiasmata = data.chiasmata;
clear data;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% select examples
chiasma = [
    19453, 1;
    158, 2];
chiasma = array2table( ...
    chiasma, 'VariableNames', {'axonId', 'chiasmaId'});

%%
paramForChiasma = transpose(horzcat( ...
    fieldnames(chiasmaParam), struct2cell(chiasmaParam)));
paramForChiasma = Util.modifyStruct(param, paramForChiasma{:});

for curIdx = 1:size(chiasma, 1)
    curAxonId = chiasma.axonId(curIdx);
    curChiasmaId = chiasma.chiasmaId(curIdx);
    
    curAxon = axons(curAxonId);
    curChiasmata = chiasmata{curAxonId};
    
    curNodeId = curChiasmata.ccCenterIdx(curChiasmaId);
    curExitNodeIds = curChiasmata.queryIdx{curChiasmaId};
    
   [~, curInnerNodeIds] = ...
        connectEM.Chiasma.restrictToSphere( ...
            param, curAxon, curNodeId, chiasmaParam.sphereRadiusInner);
    
    % add immediate neighbours as well
    curNodeIds = any(ismember(curAxon.edges, curInnerNodeIds), 2);
    curNodeIds = unique(curAxon.edges(curNodeIds, :));
    curNodeCount = numel(curNodeIds);
    
    % restrict
    curAxon.nodes = curAxon.nodes(curNodeIds, :);
   [~, curAxon.edges] = ismember(curAxon.edges, curNodeIds);
    curAxon.edges(~all(curAxon.edges, 2), :) = [];
    
   [~, curExitNodeIds] = ismember(curExitNodeIds, curNodeIds);
    curExitNodeIds = sort(reshape(curExitNodeIds, 1, 4));
    assert(all(curExitNodeIds));
    
    %%
    curIdA = curExitNodeIds(1);
    for curIdB = curExitNodeIds(2:end)
        curIdsAB = horzcat(curIdA, curIdB);
        curIdsCD = setdiff(curExitNodeIds, curIdsAB);
        
        % add loops
        curEdges = curAxon.edges;
        curEdges = cat(1, curEdges, [curIdA, curIdB; curIdsCD]);
        assert(issorted(curEdges, 2));
        
        curAdj = sparse( ...
            curEdges(:, 2), curEdges(:, 1), ...
            true, curNodeCount, curNodeCount);
        assert(graphconncomp(curAdj, 'Directed', false) == 1);

        % search for loops
        for curEdgeId = 1:size(curEdges, 1)
            curAdj = curEdges;
            curAdj(curEdgeId, :) = [];
            assert(issorted(curAdj, 2));
        
            curAdj = sparse( ...
                curAdj(:, 2), curAdj(:, 1), ...
                true, curNodeCount, curNodeCount);
           [~, curComps] = graphconncomp(curAdj, 'Directed', false);
           
            curCompAB = unique(curComps(curIdsAB));
            curCompCD = unique(curComps(curIdsCD));
            
            assert(isscalar(curCompAB));
            assert(isscalar(curCompCD));
            
            if curCompAB ~= curCompCD
                fprintf('\nFound bridge:\n');
                fprintf('(%d, %d) and (%d, %d)\n', curIdsAB, curIdsCD);
                fprintf('Edge: (%d, %d)\n', curEdges(curEdgeId, :));
            end
        end
    end
end