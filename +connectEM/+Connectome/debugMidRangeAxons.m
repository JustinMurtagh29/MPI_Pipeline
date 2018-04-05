% Check what's going on with the whacky axons that have
% * at least ten synapses
% * a spine synapse fraction between 30 % and 70%
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
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

syn = load(synFile);
conn = load(connFile);

%% Build synapse table
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

%% Limit to axons with enough synapses
axonMask = (conn.axonMeta.synCount >= minSynCount);
axonMeta = conn.axonMeta(axonMask, :);
clear axonMask;

axonMeta.spineSynFrac = ...
    axonMeta.spineSynCount ...
 ./ axonMeta.synCount;

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
numDigits = ceil(log10(1 + numel(randIds)));

mkdir(outDir);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(randIds)
    curMeta = axonMeta(randIds(curIdx), :);
    
    curAxon = conn.axons(curMeta.id);
    curSynT = synT(synT.preAggloId == curMeta.id, :);
    
    curSyns = syn.synapses(curSynT.id, :);
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
            @(id, tag) sprintf( ...
                'Synapse %d (%s)', id, spineTag{1 + tag}), ...
            curSynT.id, curSynT.isSpine, 'UniformOutput', false)];
    
    curSkel = Skeleton.fromMST(curNodes, param.raw.voxelSize, skel);
    curSkel.names = curNames;
    
    curSkelName = fullfile(outDir, sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curMeta.id));
    curSkel.write(curSkelName);
end