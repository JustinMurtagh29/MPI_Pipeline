% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_18_b_large_mergedFlightPaths.mat');

sampleNmlFile = '';
sampleNmlCount = 100;

% Path to NML file with annotated agglomerates
annNmlFile = fullfile( ...
    fileparts(mfilename('fullpath')), ...
    'annotations', 'axon-glia-classification.nml');

filterThresh = 0.3;
filterNmlFile = '/home/amotta/Desktop/glia-in-axon-candidates.nml';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segMeta = struct;
segMeta.maxSegId = Seg.Global.getMaxSegId(param);
segMeta.points = Seg.Global.getSegToPointMap(param);
segMeta = connectEM.addSegmentClassInformation(param, segMeta);

axon = load(axonFile);

%% Calculate per-axon glia score
axons = axon.axons;

gliaScore = cellfun(@(ids) median( ...
    segMeta.gliaProb(ids), 'omitnan'), axons);

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
    randIds = reshape(randIds, [], 1);
	
    numDigits = ceil(log10(1 + numel(randIds)));
    
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));
    
    for curIdx = 1:numel(randIds)
        curId = randIds(curIdx);
        curProb = randProbs(curIdx);
        curSegIds = axon.axons{curId};
        
        skel = Skeleton.fromMST( ...
            segMeta.points(curSegIds, :), ...
            param.raw.voxelSize, skel);
        skel.names{end} = sprintf( ...
            '%0*d. Axon %d (%.1f %% glia)', ...
            numDigits, curIdx, curId, 100 * curProb);
    end
    
    skel.write(sampleNmlFile);
end

%% Parse annotations
if ~isempty(annNmlFile)
    nml = slurpNml(annNmlFile);
    names = nml.things.name;
    
    aggloIds = regexpi(names, '^\d+\. Axon (\d+)', 'tokens', 'once');
    tags = regexpi(names, '\((\w+)\)$', 'tokens', 'once');
    
    % Sanity check
    assert(all(cellfun(@numel, aggloIds)));
    
    ann = table;
    ann.id = cellfun(@str2double, cat(1, aggloIds{:}));
    
    ann(cellfun(@isempty, tags), :) = [];
    ann.tag = categorical(lower(cat(1, tags{:})));
    
    ann.gliaScore = gliaScore(ann.id);
    ann.isAxon = ann.tag == 'axon';
    
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [460, 430];
    
    ax = axes(fig);
    axis(ax, 'square');
    hold(ax, 'on');
    
    binEdges = linspace(0, 1, 11);
    plotFunc = @(values, varargin) histogram( ...
        ax, values, 'BinEdges', linspace(0, 1, 11), ...
        'DisplayStyle', 'stairs', 'LineWidth', 2, varargin{:});
    plotFunc(ann.gliaScore, 'EdgeColor', 'black');
    plotFunc(ann.gliaScore(ann.isAxon));
    
    ax.TickDir = 'out';
    xlim(ax, [0, 1]);
    
    title(ax, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontWeight', 'norma', 'FontSize', 10);
end

%% Export candidates
if ~isempty(filterNmlFile)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));
    
    candIds = find(gliaScore >= filterThresh);
    numDigits = ceil(log10(1 + numel(candIds)));
    candAgglos = axons(candIds);
    
    % sort by size
   [~, sortIds] = sort(cellfun(@numel, candAgglos), 'descend');
    candAgglos = candAgglos(sortIds);
    candIds = candIds(sortIds);
    
    candAgglos = cellfun( ...
        @(ids) segMeta.points(ids, :), ...
        candAgglos, 'UniformOutput', false);
    skel = Skeleton.fromMST(candAgglos, param.raw.voxelSize, skel);
    
    skel.names = arrayfun( ...
        @(idx, id) sprintf('%0*d. Axon %d', numDigits, idx, id), ...
        transpose(1:numel(candIds)), candIds(:), 'UniformOutput', false);
    skel.write(filterNmlFile);
end
