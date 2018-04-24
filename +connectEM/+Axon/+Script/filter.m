% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');

sampleNmlFile = '/home/amotta/Desktop/median-glia.nml';
sampleNmlCount = 100;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segmentMeta = struct;
segmentMeta.maxSegId = Seg.Global.getMaxSegId(param);
segmentMeta = connectEM.addSegmentClassInformation(param, segmentMeta);

axon = load(axonFile);

%% Calculate per-axon glia score
axonIds = find(axon.indBigAxons);

axons = axon.axons(axonIds);
axons = Agglo.fromSuperAgglo(axons);

gliaScore = segmentMeta.gliaProb;
gliaScore = cellfun(@(ids) median(gliaScore(ids), 'omitnan'), axons);

%% Histogram of glia score
binEdges = linspace(0, 1, 51);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [460, 430];

ax = axes(fig);

histogram(ax, ...
    gliaScore, binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.Box = 'off';
ax.TickDir = 'out';
ax.YScale = 'log';

ax.XLim = binEdges([1, end]);
ax.YLim(1) = 0;

axis(ax, 'square');
xlabel(ax, 'Median glia probability');
ylabel(ax, 'Axons');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Export examples
if ~isempty(sampleNmlFile)
    % Select axons to export
    rng(0);
    randProbs = rand(sampleNmlCount, 1);
    
   [~, randIds] = pdist2( ...
        gliaScore, randProbs, ...
        'euclidean', 'Smallest', 1);
    randIds = axonIds(randIds(:));
	
    numDigits = ceil(log10(1 + numel(randIds)));
    
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));
    
    for curIdx = 1:numel(randIds)
        curId = randIds(curIdx);
        curProb = randProbs(curIdx);
        curAxon = axon.axons(curId);
        
        curName = sprintf( ...
            'Axon %d (%.1f %% glia)', ...
            curId, 100 * curProb);
        skel = skel.addTree( ...
            curName, curAxon.nodes(:, 1:3), curAxon.edges);
    end
    
    skel.write(sampleNmlFile);
end
