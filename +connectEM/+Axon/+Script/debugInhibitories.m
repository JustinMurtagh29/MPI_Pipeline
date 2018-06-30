% Search for shaft-preferring excitatory axons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

%% Prepare data
axonMeta = conn.axonMeta;
axonMeta.fullSpineSynFrac = ...
    axonMeta.fullSpineSynCount ...
 ./ axonMeta.fullSynCount;
axonMeta.fullPriSpineSynFrac = ...
    axonMeta.fullPriSpineSynCount ...
 ./ axonMeta.fullSynCount;

axonMeta.fullSecSpineSynFrac = ...
    axonMeta.fullSpineSynCount ...
  - axonMeta.fullPriSpineSynCount;
axonMeta.fullSecSpineSynFrac = ...
    axonMeta.fullSecSpineSynFrac ...
 ./ axonMeta.fullSynCount;

excSomaIds = find( ...
    conn.denMeta.targetClass == 'Somata' ...
  & not(conn.denMeta.isInterneuron));

excSomaConn = conn.connectome;
excSomaConn(~ismember(excSomaConn.edges(:, 2), excSomaIds), :) = [];
excSomaConn.synCount = cellfun(@numel, excSomaConn.synIdx);

axonMeta.excSomaSynCount = accumarray( ...
    excSomaConn.edges(:, 1), ...
    excSomaConn.synCount, ...
   [height(conn.axonMeta), 1]);
axonMeta.excSomaSynFrac = ...
    axonMeta.excSomaSynCount ...
 ./ axonMeta.synCount;

axonMeta(axonMeta.synCount < minSynPre, :) = [];

%% Prepare for plot
axonMeta.exclusion = zeros(height(axonMeta), 1);
axonMeta.exclusion(axonMeta.fullSpineSynFrac > 0.5) = 1;
axonMeta.exclusion(axonMeta.fullPriSpineSynFrac > 0.5) = 2;
axonMeta.exclusion(axonMeta.fullSecSpineSynFrac > 0.3) = 3;
axonMeta.exclusion(axonMeta.excSomaSynFrac > 0.1) = 4;

colors = get(groot, 'defaultAxesColorOrder');
axonMeta.color = colors(1 + axonMeta.exclusion, :);

%% Plot
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [675, 330];

ax = subplot(1, 3, 1);
scatter(ax, ...
    axonMeta.fullPriSpineSynFrac, ...
    axonMeta.excSomaSynFrac, ...
    36, axonMeta.color, '.');
xlabel(ax, { ...
    'Fraction of synapses being'; ...
    'primary spine innervations'});
ylabel(ax, { ...
    'Fraction of synapses'; ...
    'onto excitatory somata'});

ax = subplot(1, 3, 2);
scatter(ax, ...
    axonMeta.fullPriSpineSynFrac, ...
    axonMeta.fullSecSpineSynFrac, ...
    36, axonMeta.color, '.');
ylabel(ax, { ...
    'Fraction of synapses being'; ...
    'secondary spine innervations'});

ax = subplot(1, 3, 3);
scatter(ax, ...
    axonMeta.fullPriSpineSynFrac, ...
    axonMeta.fullSpineSynFrac, ...
    36, axonMeta.color, '.');
ylabel(ax, { ...
    'Fraction of synapses'; ...
    'onto spines'});

set(fig.Children, ...
    'XLim', [0, 1], ...
    'YLim', [0, 1], ...
    'TickDir', 'out', ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Search for potential shaft-preferring excitatory axons
rng(0);
axonIds = axonMeta.id(~axonMeta.exclusion);
axonIds = randperm(numel(excCandIds));
axonIds = excCandIds(axonIds);
axonIds = axonIds(1:25);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

numDigits = ceil(log10(1 + numel(axonIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(axonIds)
    curId = conn.axonMeta.id(axonIds(curIdx));
    curSynT = synT(synT.preAggloId == curId, :);
    
    curAxon = conn.axons(curId);
    curSynapses = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynT.id), ...
        syn.synapses.postsynId(curSynT.id), ...
        'UniformOutput', false);
    
    curNodes = [curAxon; curSynapses];
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id, type) sprintf('Synapse %d (%s)', id, type), ...
        curSynT.id, curSynT.type, 'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
