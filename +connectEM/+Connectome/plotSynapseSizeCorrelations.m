% TODO(amotta)
% * `contactArea` in `connectomeMeta` ranges up to ~120. That's true even
%   for spine-head filtered synapses. Find out what's going on: Do I
%   misunderstand the data structure or the physical units? Is there a bug?
%   Is synapse agglomeration off? Soma segments classified as spine?
%
% Written by`
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_with_den_meta.mat');

% set this variable to debug
debugDir = '';

info = Util.runInfo();

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
conn = load(connFile);
syn = load(synFile);

% for debugging
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

%% limit synapses
synT = table;
synT.id = cell2mat(conn.connectome.synIdx);
synT.area = cell2mat(conn.connectomeMeta.contactArea);
synT.isSpine = syn.isSpineSyn(synT.id);

synT.preAggloId = repelem( ...
    conn.connectome.edges(:, 1), ...
    cellfun(@numel, conn.connectome.synIdx));

synT.postAggloId = repelem( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx));

% limit to spine synapses
synT(~synT.isSpine, :) = [];
synT.isSpine = [];

% remove duplicate entries
[~, uniRows, uniCount] = unique(synT.id);
synT = synT(uniRows, :);

% remove synapses occuring multiple times
% (i.e., between at least two different pairs of neurites)
synT.occurences = accumarray(uniCount, 1);
synT(synT.occurences > 1, :) = [];
synT.occurences = [];

% remove synapses whose size is obviously wrong
synT(synT.area > 1.5, :) = [];

%% look at doubly coupled neurites
[dupNeurites, ~, uniRows] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

dupNeurites.synIds = accumarray( ...
    uniRows, synT.id, [], @(synIds) {transpose(synIds)});

uniCount = accumarray(uniRows, 1);
dupNeurites(uniCount ~= 2, :) = [];

% we now know that there are exactly two synapses
dupNeurites.synIds = cell2mat(dupNeurites.synIds);

% load synapse areas
[~, synRows] = ismember(dupNeurites.synIds, synT.id);
dupNeurites.synAreas = synT.area(synRows);

flipMask = ...
    dupNeurites.synAreas(:, 2) ...
  < dupNeurites.synAreas(:, 1);

dupNeurites.synIds(flipMask, :) = ...
    fliplr(dupNeurites.synIds(flipMask, :));
dupNeurites.synAreas(flipMask, :) = ...
    fliplr(dupNeurites.synAreas(flipMask, :));

%%
fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    dupNeurites.synAreas(:, 2), ...
    dupNeurites.synAreas(:, 1), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

% do the fit
% taken from `matlab/+Analysis/+Script/bartolEtAl2015eLife.m` in `amotta`
xLog = log10(dupNeurites.synAreas(:, 2));
yLog = log10(dupNeurites.synAreas(:, 1));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('y = %.2f x^{%.2f}', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));
plot(fitRange, fitRange, 'k--');

title({'Same-axon same-dendrite spine synapses'; info.git_repos{1}.hash});
legend(rawName, fitName, 'Location', 'NorthWest');

%% look at highly inconsistent pairs
% search for pairs which differ by more than a factor of 10
randDupNeurites = find( ...
    dupNeurites.synAreas(:, 2) ...
 ./ dupNeurites.synAreas(:, 1) > 10);

rng(0);
randDupNeurites = randDupNeurites( ...
    randperm(numel(randDupNeurites), 25));
randDupNeurites = dupNeurites(randDupNeurites, :);

for curIdx = 1:size(randDupNeurites, 1)
    curRow = randDupNeurites(curIdx, :);
    
    curSynAgglos = cellfun(@vertcat, ...
        syn.synapses.presynId(curRow.synIds), ....
        syn.synapses.postsynId(curRow.synIds), ...
        'UniformOutput', false);
    
    curAgglos = [ ...
        conn.axons(curRow.preAggloId); ...
        conn.dendrites(curRow.postAggloId); ...
        curSynAgglos];
    curNodes = cellfun( ...
        @(segIds) points(segIds, :), ...
        curAgglos, 'UniformOutput', false);
    
    curNames = { ...
        sprintf('Axon #%d', curRow.preAggloId); ...
        sprintf('Dendrite #%d', curRow.postAggloId)};
    curNames = [curNames; arrayfun(@(synId) ....
        sprintf('Synapse #%d', synId), ...
        curRow.synIds(:), 'UniformOutput', false)]; %#ok
    
    curSkel = Skeleton.fromMST(curNodes, param.raw.voxelSize);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel.names = curNames;
    
    curSkelFile = sprintf( ...
        '%d_axon-%d_dendrite_%d.nml', ...
        curIdx, curRow.preAggloId, curRow.postAggloId);
    curSkelFile = fullfile(debugDir, curSkelFile);
    
    if isempty(debugDir); continue; end
    curSkel.write(curSkelFile);
end

%% plot distribution of synapse size
% plot distribution
fig = figure();
ax = axes(fig);

histogram(ax, synT.area, linspace(0, 1.5, 51));
xlabel(ax, 'Axon-spine interface (µm²)');
ylabel(ax, 'Spine synapses');
ax.TickDir = 'out';

annotation(...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Spine synapse size distribution';
        info.git_repos{1}.hash});
    
fig.Position(3:4) = [820, 475];
