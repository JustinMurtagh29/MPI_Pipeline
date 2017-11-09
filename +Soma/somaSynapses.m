% script to exclude soma synapses on the exits of the soma
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();
doPlot = false; %#ok<*UNRCH>


%% load data
 
p = Gaba.getSegParameters('ex145_ROI2017');

% graph
[graph, segmentMeta, borderMeta] = Seg.IO.loadGraph(p, false);

% somas
m = load(p.agglo.somaFile);
somaAgglos = m.somaAgglos(:,1);
numSomas = length(somaAgglos);

% synapses
m = load(p.connectome.synFile);
synapses = m.synapses;
boutons = m.boutons;
isSpineSyn = m.isSpineSyn;


%% get soma synapses

maxSomaId = max(cellfun(@max, somaAgglos(:,1)));
postsynSomaLUT = L4.Agglo.buildLUT(synapses.postsynId, maxSomaId);
somaSynIdxAll = cellfun(@(x)setdiff(postsynSomaLUT(x), 0), ...
    somaAgglos(:,1), 'uni', 0);
somaSynIdx = somaSynIdxAll;

% soma syn coms
somaSynComs = cell(numSomas, 1);
for i = 1:length(somaSynIdx)
    somaSynComs{i} = round(cell2mat( ...
        cellfun(@(x)mean(borderMeta.borderCoM(x, :), 1), ...
        synapses.edgeIdx(somaSynIdx{i}), 'uni', 0)));
end

synScoresT.synT = -1;
synScoresT.otherSynT = 1.5;
toDiscard1 = cell(numSomas, 1);
for i = 1:numSomas
    [~, toDiscard1{i}] = L4.Soma.somaSynapsePostProcessing( ...
        synapses, boutons, somaSynIdx{i}, ...
        isSpineSyn, graph.synScores, synScoresT );
end
% discardedSynIdx = cellfun(@(x,y)x(y), somaSynIdx, toDiscard, 'uni', 0);
discardedSynComs = cellfun(@(x, idx) x(idx, :), somaSynComs, toDiscard1, ...
    'uni', 0);
somaSynIdx = cellfun(@(x,y)x(~y), somaSynIdx, toDiscard1, 'uni', 0);
somaSynComs = cellfun(@(x, idx) x(~idx, :), somaSynComs, toDiscard1, ...
    'uni', 0);

% somaSynEdgeIdx = cellfun(@(x)cell2mat(synapses.edgeIdx(x)), ...
%     somaSynIdx, 'uni', 0);

[~, ~, toDiscard2] = Soma.morphPostProc(p, somaAgglos, ...
    segmentMeta.point, somaSynComs);

% discardedSynIdx2 = cellfun(@(x,y)x(y), somaSynIdx, toDiscardSyn, 'uni', 0);
discardedSynComs2 = cellfun(@(x, idx) x(idx, :), somaSynComs, toDiscard2, ...
    'uni', 0);

somaSynIdx = cellfun(@(x,y)x(~y), somaSynIdx, toDiscard2, 'uni', 0);
somaSynComs = cellfun(@(x, idx) x(~idx, :), somaSynComs, toDiscard2, ...
    'uni', 0);


%% save results

[outFolder, outFile] = fileparts(p.agglo.somaFile);
outFile = fullfile(outFolder, [outFile '_synapses.mat']);
if ~exist(outFile, 'file')
    save(outFile, 'info', 'somaSynIdxAll', 'toDiscard1', 'toDiscard2', ...
        'somaSynIdx', 'somaSynComs');
else
    Util.log('Output file %s already exists. Save the file manually.', ...
        outFile);
end


%% some visualization

if doPlot
    
    % spine synapse exclusion
    idx = cellfun(@isempty, discardedSynComs);
    idx = find(~idx);
    skel = Skeleton.fromMST(discardedSynComs(idx));
    for i = 1:length(idx)
        skel.names{i} = sprintf('Soma%02d', idx(i));
    end
    skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');
    skel.write('SomaSpineSynExclusion.nml');
    
    % exit synapse exclusion
    idx = cellfun(@isempty, discardedSynComs2);
    idx = find(~idx);
    skel = Skeleton.fromMST(discardedSynComs2(idx));
    for i = 1:length(idx)
        skel.names{i} = sprintf('Soma%02d', idx(i));
    end
    skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');
    skel.write('SomaExitSynExclusion.nml');
    
    % soma synapses
    idx = cellfun(@isempty, somaSynComs);
    idx = find(~idx);
    skel = Skeleton.fromMST(somaSynComs(idx));
    for i = 1:length(idx)
        skel.names{i} = sprintf('Soma%02d', idx(i));
    end
    skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');
    skel.write('SomaSynComs.nml');
end

