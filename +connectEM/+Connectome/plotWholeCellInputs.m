% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/tmpscratch/amotta/l4/2018-02-23-whole-cell-input-distributions/20180223T132250-run';

synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segMass = Seg.Global.getSegToSizeMap(param);
segCentroids = Seg.Global.getSegToCentroidMap(param);

syn = load(synFile);
dend = load(dendFile);

[~, connName] = fileparts(connFile);
conn = connectEM.Connectome.load(param, connName);

%% establish correspondence between connectome and whole cells
somaAggloIds = find(conn.denMeta.targetClass == 'Somata');
somaAgglos = conn.dendrites(somaAggloIds);
somaLUT = Agglo.buildLUT(maxSegId, somaAgglos);

dendAggloIds = find(conn.denMeta.targetClass == 'WholeCell');
dendAgglos = conn.dendrites(dendAggloIds);
dendLUT = Agglo.buildLUT(maxSegId, dendAgglos);

wcT = table;
wcT.aggloId = dend.indWholeCells;
wcT.segIds = dend.dendAgglos(wcT.aggloId);

% find soma
wcT.somaId = cellfun( ...
    @(segIds) setdiff(somaLUT(segIds), 0), ...
    wcT.segIds, 'UniformOutput', false);

wcT(~cellfun(@isscalar, wcT.somaId), :) = [];
wcT.somaId = cell2mat(wcT.somaId);

% find dendrite(s)
wcT.dendId = cellfun( ...
    @(segIds) setdiff(dendLUT(segIds), 0), ...
    wcT.segIds, 'UniformOutput', false);

wcT(~cellfun(@isscalar, wcT.dendId), :) = [];
wcT.dendId = cell2mat(wcT.dendId);

% add soma to segment list
wcT.segIds = cellfun( ...
    @(segIds, somaSegIds) unique(cat(1, segIds, somaSegIds)), ...
    wcT.segIds, somaAgglos(wcT.somaId), 'UniformOutput', false);

%% calculate distance to soma surface
wcT.nodeDists = cell(size(wcT.segIds));

for curIdx = 1:size(wcT, 1)
    curSegIds = wcT.segIds{curIdx};
    
    curSomaSegIds = somaAgglos{wcT.somaId(curIdx)};
   [~, curSomaSegIds] = ismember(curSomaSegIds, curSegIds);
    assert(all(curSomaSegIds));
    
    % build weights
    curDists = segCentroids(curSegIds, :);
    curDists = curDists .* param.raw.voxelSize;
    curDists = pdist(curDists);
    assert(all(curDists));
    
    curDists = squareform(curDists);
    
    % NOTE(amotta): Set distance between somatic segments to zero.
    %
    % True zeros cannot be used (not even with the additional 'Weights'
    % option of `graphminspantree`) since the output matrix will contain
    % zeros for both absent and somatic edges.
    curDists(curSomaSegIds, curSomaSegIds) = eps;
    
    curAdj = graphminspantree( ...
        sparse(curDists), curSomaSegIds(1));
    curAdj = logical(curAdj);
    
    % now we can use true zeros
    curDists(curSomaSegIds, curSomaSegIds) = 0;

    curDists = graphshortestpath( ...
        curAdj, curSomaSegIds(1), ...
        'Weights', curDists(curAdj), ...
        'Directed', false);
    
    % sanity checks
    curSomaMask = false(size(curSegIds));
    curSomaMask(curSomaSegIds) = true;
    
    assert(~any(curDists(curSomaMask)));
    assert(all(curDists(~curSomaMask)));
    
    curDists = reshape(curDists, [], 1);
    wcT.nodeDists{curIdx} = curDists;
end

%% collect input synapses
wcT.synapses = cell(size(wcT.aggloId));

for curIdx = 1:size(wcT, 1)
    curPostIds = [ ...
        somaAggloIds(wcT.somaId(curIdx)), ...
        dendAggloIds(wcT.dendId(curIdx))];
    
    curConnRows = ismember( ...
        conn.connectome.edges(:, 2), curPostIds);
    curConnRows = conn.connectome(curConnRows, :);
    
    curSynT = table;
    curSynT.id = cell2mat(curConnRows.synIdx);
    
    curSynT.segIdx = syn.synapses.postsynId(curSynT.id);
    curSynT.segIdx = cellfun( ...
        @(ids) intersect(ids, wcT.segIds{curIdx}), ...
        curSynT.segIdx, 'UniformOutput', false);
    
    % assign synapse to largest segment
   [~, curSegIdx] = cellfun(@(ids) max(segMass(ids)), curSynT.segIdx);
    curSegIds = arrayfun(@(ids, idx) ids{1}(idx), curSynT.segIdx, curSegIdx);
    
    % translate to relative 
   [~, curSynT.segIdx] = ismember(double(curSegIds), wcT.segIds{curIdx});
    
    curSynT.axonId = repelem( ...
        curConnRows.edges(:, 1), ...
        cellfun(@numel, curConnRows.synIdx));
    
    wcT.synapses{curIdx} = curSynT;
end

%% split axons into exc. and inh.
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
 ./ conn.axonMeta.synCount;

conn.axonMeta.isExc = (conn.axonMeta.spineSynFrac > 0.5);
conn.axonMeta.isInh = ~conn.axonMeta.isExc;

%% plotting
for curIdx = 1:size(wcT, 1)
    curFig = figure('visible', 'off');
    curFig.Position(3:4) = [860, 480];
    
    curSyns = wcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.isSoma = logical(somaLUT( ...
        wcT.segIds{curIdx}(curSyns.segIdx)));
    
    curSyns.isExc = conn.axonMeta.isExc(curSyns.axonId);
    curSyns.isInh = conn.axonMeta.isInh(curSyns.axonId);
    curSyns.isTc = conn.axonMeta.isThalamocortical(curSyns.axonId);
    
    curSyns.dist = wcT.nodeDists{curIdx}(curSyns.segIdx);
    curSyns.dist = curSyns.dist / 1E3;
    
    % move soma synapses to separate bin
    curSyns.dist(curSyns.isSoma) = -eps;
    
    curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
    curBinEdges = -10:10:curMaxDist;
    
    curPlot = @(ax, data) ...
        histogram( ...
            ax, data, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
    curAx = subplot(2, 1, 1);
    hold(curAx, 'on');
    
	curPlot(curAx, curSyns.dist);
    curPlot(curAx, curSyns.dist(curSyns.isSpine));
    curPlot(curAx, curSyns.dist(curSyns.isSoma));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Synapses');
    
    title(curAx, { ...
        sprintf('Inputs onto whole cell %d', curIdx); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    curLeg = legend(curAx, ...
        'All', 'Onto spines', 'Onto soma', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    
    curAx = subplot(2, 1, 2);
    hold(curAx, 'on');
    
    curPlot(curAx, curSyns.dist(curSyns.isExc));
    curPlot(curAx, curSyns.dist(curSyns.isInh));
    curPlot(curAx, curSyns.dist(curSyns.isTc));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlabel(curAx, 'Distance to soma (µm)');
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Synapses');
    
    curLeg = legend(curAx, ...
        'Excitatory', ...
        'Inhibitory', ...
        'Thalamocortical', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    
    curFigName = sprintf('input-distribution_whole-cell-%d.png', curIdx);
    curFigName = fullfile(outputDir, curFigName);
    
    export_fig('-r172', curFigName, curFig);
    clear curFig;
end