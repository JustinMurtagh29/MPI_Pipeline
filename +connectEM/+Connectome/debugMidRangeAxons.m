% Check what's going on with the whacky axons that have
% * at least ten synapses
% * a spine synapse fraction between 30 % and 70%
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
outDir = '/home/amotta/Desktop/mid-range-axons';

minSynCount = 10;
minSpineFrac = 0.3;
maxSpineFrac = 0.7;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

syn = load(synFile);
conn = load(connFile);

%% Collect synapses
shLUT = (Agglo.buildLUT(maxSegId, shAgglos) ~= 0);

syn.synapses.id = reshape( ...
    1:size(syn.synapses, 1), [], 1);
syn.synapses.ontoSpine = cellfun( ...
    @(segIds) any(shLUT(segIds)), ...
    syn.synapses.postsynId);
clear shLUT;

% This time we count a synapse even when none of the postsynaptic segments
% is in a dendrite agglomerates. This should give a more accurate picture
% of the axonal outputs.
axonLUT = Agglo.buildLUT(maxSegId, conn.axons);

synapses = syn.synapses;
synapses.axonId = cellfun( ...
    @(segIds) setdiff(axonLUT(segIds), 0), ...
    synapses.presynId, 'UniformOutput', false);
clear axonLUT;

% Remove synapses which have
% * no presynaptic axon at all
% * multiple axons on presynaptic side (?!)
synapses(~cellfun(@isscalar, synapses.axonId), :) = [];
synapses.axonId = cell2mat(synapses.axonId);

axonMeta = conn.axonMeta;
axonMeta.fullSynCount = accumarray( ...
    synapses.axonId, 1, size(axonMeta.id));
axonMeta.fullSpineSynCount = accumarray( ...
    synapses.axonId, synapses.ontoSpine, size(axonMeta.id));

axonMeta.synIds = accumarray( ...
    conn.connectome.edges(:, 1), ...
    transpose(1:size(conn.connectome, 1)), size(axonMeta.id), ...
    @(r) {cell2mat(conn.connectome.synIdx(r))}, {zeros(0, 1)});
axonMeta.fullSynIds = accumarray( ...
    synapses.axonId, synapses.id, size(axonMeta.id), ...
    @(synIds) {synIds}, {zeros(0, 1)});

% Sanity checks
assert(all(cellfun(@numel, axonMeta.synIds) == axonMeta.synCount));
assert(all(cellfun(@numel, axonMeta.fullSynIds) == axonMeta.fullSynCount));

axonMeta.spineSynFrac = ...
    axonMeta.spineSynCount ...
 ./ axonMeta.synCount;
axonMeta.fullSpineSynFrac = ...
    axonMeta.fullSpineSynCount ...
 ./ axonMeta.fullSynCount;

clear synapses;

%% Plot spine synapse fractions
binEdges = linspace(0, 1, 21);
axonMask = (axonMeta.synCount >= minSynCount);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square')
hold(ax, 'on');

histogram(ax, ...
    axonMeta.spineSynFrac(axonMask), binEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(ax, ...
    axonMeta.fullSpineSynFrac(axonMask), binEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));
xlabel(ax, 'Spine synapse fraction');
ylabel(ax, 'Axons');

legend(ax, ...
    'Synapses onto dendrite agglos', ...
    'All synapses', 'Location', 'NorthWest');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot distribution of spine synapse fraction
fig = figure;
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

histogram(ax, ...
    axonMeta.spineSynFrac, ...
    linspace(0, 1, 21), ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

%% Pick random examples (corrected for frequency)
axonMeta(axonMeta.spineSynFrac < minSpineFrac, :) = [];
axonMeta(axonMeta.spineSynFrac > maxSpineFrac, :) = [];

rng(0);
axonMeta = sortrows(axonMeta, 'id');
axonMeta = axonMeta(randperm(size(axonMeta, 1)), :);

rng(0);
randIds = rand(50, 1);
randIds = randIds .* (maxSpineFrac - minSpineFrac);
randIds = randIds + minSpineFrac;

randIds = pdist2(axonMeta.spineSynFrac, randIds);
[~, randIds] = min(randIds, [], 1);

%% Export to webKNOSSOS
spineTag = {'shaft', 'spine'};
connTag = {'not in conn.', 'in conn.'};
numDigits = ceil(log10(1 + numel(randIds)));

mkdir(outDir);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(randIds)
    curMeta = axonMeta(randIds(curIdx), :);
    
    curAxon = conn.axons(curMeta.id);
    curSynIds = curMeta.fullSynIds{1};
    curSynOntoSpine = syn.synapses.ontoSpine(curSynIds);
    curSynInConn = ismember(curSynIds, curMeta.synIds{1});
    
    curSyns = syn.synapses(curSynIds, :);
    curSyns = cellfun(@vertcat, ...
        curSyns.presynId, curSyns.postsynId, ...
        'UniformOutput', false);
    
    curNodes = vertcat(curAxon, curSyns);
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curNames = [ ...
       {sprintf( ...
            'Axon %d (%.1f %% of synapses onto spines)', ...
            curMeta.id, 100 * curMeta.spineSynFrac)}; ...
        arrayfun( ...
            @(id, spine, conn) sprintf( ...
                'Synapse %d (%s, %s)', ...
                id, spineTag{1 + spine}, connTag{1 + conn}), ...
            curSynIds, curSynOntoSpine, curSynInConn, ...
            'UniformOutput', false)];
    
    curSkel = Skeleton.fromMST(curNodes, param.raw.voxelSize, skel);
    curSkel.names = curNames;
    
    curSkelName = fullfile(outDir, sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curMeta.id));
    curSkel.write(curSkelName);
end