% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

[interSynDir, interSynFile] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse.mat', interSynFile);
interSynFile = fullfile(interSynDir, interSynFile);
clear interSynDir;

% set this variable to debug
debugDir = '';

info = Util.runInfo();

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
conn = load(connFile);
syn = load(synFile);
interSyn = load(interSynFile);

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
% (i.e., between two or more pairs of neurites)
synT.occurences = accumarray(uniCount, 1);
synT(synT.occurences > 1, :) = [];
synT.occurences = [];

%% plot distribution of synapse size
% plot distribution
fig = figure();
ax = axes(fig);

histogram(ax, synT.area, 51);
xlabel(ax, 'Axon-spine interface (µm²)');
ylabel(ax, 'Spine synapses');
ax.TickDir = 'out';

title( ...
   {'Spine synapse size distribution'; info.git_repos{1}.hash}, ...
    'FontWeight', 'norma', 'FontSize', 10);
    
fig.Position(3:4) = [820, 475];

%% plot histogram over no. of synapse per neurite pair
% precompute
[~, ~, neuriteCoupling] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');
neuriteCoupling = accumarray(neuriteCoupling, 1);

% pot
fig = figure();
ax = axes(fig);

histogram(ax, neuriteCoupling, 'LineWidth', 2, 'DisplayStyle', 'stairs');
title(ax, info.git_repos{1}.hash, 'FontWeight', 'normal', 'FontSize', 10);
xlabel(ax, 'Spine synapses per connection');
ylabel(ax, 'Connections');

ax.YScale = 'log';
ax.TickDir = 'out';

ax.XLim = 0.5 + [0, max(neuriteCoupling)];
ax.YLim(1) = 10 ^ (-0.1);
ax.YTickLabel = arrayfun(@num2str, ax.YTick, 'UniformOutput', false);

%% ASI areas vs. degree of coupling
[~, ~, neuriteCoupling] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

neuriteSynAreas = accumarray( ...
    neuriteCoupling, synT.area, [], @(areas) {areas});
neuriteCoupling = accumarray(neuriteCoupling, 1);

neuriteCoupling = repelem( ...
    neuriteCoupling, cellfun(@numel, neuriteSynAreas));
neuriteSynAreas = cell2mat(neuriteSynAreas);

% plot
fig = figure();
ax = axes(fig);

boxplot(ax, neuriteSynAreas, neuriteCoupling);
title(ax, info.git_repos{1}.hash, 'FontWeight', 'normal', 'FontSize', 10);
xlabel(ax, 'Spine synapses per connection');
ylabel(ax, 'Axon-spine interface area (µm²)');

ax.YLim(1) = 0;
ax.TickDir = 'out';

%% ASI area variability
[~, ~, neuriteCoupling] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');
neuriteCv = accumarray( ...
    neuriteCoupling, synT.area, [], ...
    @(areas) std(areas) / mean(areas));
neuriteCoupling = accumarray(neuriteCoupling, 1);

% remove single-spine case
neuriteCv(neuriteCoupling < 2) = [];
neuriteCoupling(neuriteCoupling < 2) = [];

% plot
fig = figure();
ax = axes(fig);

boxplot(ax, neuriteCv, neuriteCoupling);
title(ax, info.git_repos{1}.hash, 'FontWeight', 'normal', 'FontSize', 10);
xlabel(ax, 'Spine synapses per connection');
ylabel(ax, 'Variability of all ASI areas (CV)');

ax.YLim(1) = 0;
ax.TickDir = 'out';

%% ASI area variability for largest two synapses
[~, ~, neuriteCoupling] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

upToTwo = @(i) i(1:min(2, numel(i)));
cvOf = @(i) std(i) / mean(i);

neuriteCv = accumarray( ...
    neuriteCoupling, synT.area, [], ...
    @(a) cvOf(upToTwo(sort(a, 'descend'))));
neuriteCoupling = accumarray(neuriteCoupling, 1);

% remove single-spine case
neuriteCv(neuriteCoupling < 2) = [];
neuriteCoupling(neuriteCoupling < 2) = [];

% plot
fig = figure();
ax = axes(fig);

boxplot(ax, neuriteCv, neuriteCoupling);
title(ax, info.git_repos{1}.hash, 'FontWeight', 'normal', 'FontSize', 10);
xlabel(ax, 'Spine synapses per connection');
ylabel(ax, 'Variability of largest two ASI areas (CV)');

ax.YLim(1) = 0;
ax.TickDir = 'out';

%% calculate baseline slope
synT = sortrows(synT, 'area', 'ascend');

rng(0);
randIds = floor(size(synT, 1) / 2);
randIds = randperm(size(synT, 1), 2 * randIds);

randIds = reshape(randIds, [], 2);
randSynAreas = synT.area(randIds);
randSynAreas = sort(randSynAreas, 2);

xLog = log10(randSynAreas(:, 2));
yLog = log10(randSynAreas(:, 1));

bl = [ones(numel(xLog), 1), xLog] \ yLog;
bl(1) = 10 ^ bl(1);

blFitF = @(x) bl(1) .* (x .^ bl(2));
blFitName = sprintf('Random pairs (y = %.2f x^{%.2f})', bl(1), bl(2));

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

[~, dupNeurites.interSynDist] = ismember( ...
    dupNeurites.preAggloId, interSyn.axonIds);
for curIdx = 1:size(dupNeurites, 1)
    curAxonIdx = dupNeurites.interSynDist(curIdx);
    curInterSynDists = interSyn.synToSynDists{curAxonIdx};
    
   [~, curSynIdx] = ismember( ...
        dupNeurites.synIds(curIdx, :), ...
        interSyn.synIds{curAxonIdx});
    dupNeurites.interSynDist(curIdx) = ...
        curInterSynDists(curSynIdx(1), curSynIdx(2));
end

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
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));
plot(fitRange, blFitF(fitRange));
plot(fitRange, fitRange, 'k--');

title( ...
   {'Same-axon same-dendrite spine synapses'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% same-axon different-dendrite pairs
rng(0);

% get rid of any synapse size preference
saddT = synT(randperm(size(synT, 1)), :);

% TODO(amotta): If axon A makes multiple synapses onto dendrite D, we
% forget about all but one synapses for the following analysis. This makes
% the code easier. But I'm not yet sure whether this introduces some kind
% of bias...
[~, uniRows] = unique(saddT(:, {'preAggloId', 'postAggloId'}), 'rows');
saddT = saddT(uniRows, :);

% find axons that occur at least twice
[dendDupIds, ~, dendDupCount] = unique(saddT.preAggloId);
dendDupCount = accumarray(dendDupCount, 1);

dendDupIds(dendDupCount < 2) = [];
dendDupCount(dendDupCount < 2) = [];

assert(issorted(saddT.preAggloId));
[~, uniRows] = ismember(dendDupIds, saddT.preAggloId);

uniRows = arrayfun( ...
    @(r, c) r + (1:(2 * floor(c / 2)))' - 1, ...
	uniRows, dendDupCount, 'UniformOutput', false);
uniRows = cell2mat(uniRows);

saddT = saddT(uniRows, :);
saddT.pairId = ceil((1:size(saddT, 1)) / 2)';

% now that we've chosen the synapses, we can sort by size
saddT = sortrows(saddT, {'preAggloId', 'pairId', 'area'});

% sanity checks
assert(all(saddT.preAggloId(1:2:end) == saddT.preAggloId(2:2:end)));
assert(all(saddT.postAggloId(1:2:end) ~= saddT.postAggloId(2:2:end)));
assert(all(saddT.area(1:2:end) <= saddT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    saddT.area(2:2:end), ...
    saddT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(saddT.area(2:2:end));
yLog = log10(saddT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Same-axon different-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% different-axon same-dendrite pairs
rng(0);

% get rid of any synapse size preference
dasdT = synT(randperm(size(synT, 1)), :);

% TODO(amotta): If axon A makes multiple synapses onto dendrite D, we
% forget about all but one synapses for the following analysis. This makes
% the code easier. But I'm not yet sure whether this introduces some kind
% of bias...
[~, uniRows] = unique(dasdT(:, {'postAggloId', 'preAggloId'}), 'rows');
dasdT = dasdT(uniRows, :);

% find axons that occur at least twice
[dendDupIds, ~, dendDupCount] = unique(dasdT.postAggloId);
dendDupCount = accumarray(dendDupCount, 1);

dendDupIds(dendDupCount < 2) = [];
dendDupCount(dendDupCount < 2) = [];

assert(issorted(dasdT.postAggloId));
[~, uniRows] = ismember(dendDupIds, dasdT.postAggloId);

uniRows = arrayfun( ...
    @(r, c) r + (1:(2 * floor(c / 2)))' - 1, ...
	uniRows, dendDupCount, 'UniformOutput', false);
uniRows = cell2mat(uniRows);

dasdT = dasdT(uniRows, :);
dasdT.pairId = ceil((1:size(dasdT, 1)) / 2)';

% now that we've chosen the synapses, we can sort by size
dasdT = sortrows(dasdT, {'postAggloId', 'pairId', 'area'});

% sanity checks
assert(all(dasdT.postAggloId(1:2:end) == dasdT.postAggloId(2:2:end)));
assert(all(dasdT.preAggloId(1:2:end) ~= dasdT.preAggloId(2:2:end)));
assert(all(dasdT.area(1:2:end) <= dasdT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    dasdT.area(2:2:end), ...
    dasdT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(dasdT.area(2:2:end));
yLog = log10(dasdT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Different-axon same-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% different-axon different-dendrite
rng(0);

% get rid of any synapse size preference
daddT = synT(randperm(size(synT, 1)), :);

[~, uniRows] = unique(daddT(:, {'preAggloId', 'postAggloId'}), 'rows');
daddT = daddT(uniRows, :);

curPos = 0;
curSynsOpen = true(size(daddT, 1), 1);
curIds = zeros(0, 2);

while curPos < height(daddT)
    curPos = curPos + 1;
    
    if ~curSynsOpen(curPos)
        continue;
    end
    
    curSyn = daddT(curPos, :);
    curPossIds = all(bsxfun(@ne, ...
        [daddT.preAggloId, daddT.postAggloId], ...
        [curSyn.preAggloId, curSyn.postAggloId]), 2);
    curPossIds = find(curPossIds & curSynsOpen);
    
    if ~isempty(curPossIds)
        curPossId = curPossIds(randi(numel(curPossIds)));
        curSynsOpen(cat(2, curPossId, curPos)) = false;
        curIds(end + 1, :) = [curPos, curPossId]; %#ok
    end
end

curIds = reshape(transpose(curIds), [], 1);
daddT = daddT(curIds, :);

daddT.pairId = ceil((1:size(daddT, 1)) / 2)';
daddT = sortrows(daddT, {'pairId', 'area'});

% sanity checks
assert(all(daddT.postAggloId(1:2:end) ~= daddT.postAggloId(2:2:end)));
assert(all(daddT.preAggloId(1:2:end) ~= daddT.preAggloId(2:2:end)));
assert(all(daddT.area(1:2:end) <= daddT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    daddT.area(2:2:end), ...
    daddT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(daddT.area(2:2:end));
yLog = log10(daddT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Different-axon different-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%%
maxCv = 1.5;

cvData = struct;
cvData(1).name = 'Random pairs';
cvData(1).data = randSynAreas;

cvData(2).name = 'Same-axon same-dendrite';
cvData(2).data = dupNeurites.synAreas;

cvData(3).name = 'Same-axon different-dendrite';
cvData(3).data = reshape(saddT.area, 2, [])';

cvData(4).name = 'Different-axon same-dendrite';
cvData(4).data = reshape(dasdT.area, 2, [])';

cvData(5).name = 'Different-axon different-dendrite';
cvData(5).data = reshape(daddT.area, 2, [])';

fig = figure();
ax = axes(fig);
hold(ax, 'on');

for curIdx = 1:numel(cvData)
    cvData(curIdx).cv = ...
        std(cvData(curIdx).data, 0, 2) ...
        ./ mean(cvData(curIdx).data, 2);
    
    cvData(curIdx).median = ...
        median(cvData(curIdx).cv);
    
   [~, cvData(curIdx).pValue] = ...
        kstest2(cvData(curIdx).cv, cvData(1).cv);
    
    if curIdx == 1; continue; end
    
    cvData(curIdx).name = sprintf( ...
        '%s (p = %.2g)', ...
        cvData(curIdx).name, ...
        cvData(curIdx).pValue);
end

for curIdx = 1:numel(cvData)
    histogram( ...
        ax, cvData(curIdx).cv, linspace(0, maxCv, 51), ...
        'Normalization', 'probability', 'DisplayStyle', 'stairs');
end

title( ...
   {'Synapse size variability'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(ax, cvData.name, 'Location', 'SouthWest');

ax.TickDir = 'out';
ax.XLim = [0, maxCv];
xlabel(ax, 'Coefficient of variation');
ylabel(ax, 'Probability');

%% estimate circuit learnedness
% TODO(amotta): justify + automatically derive
threshCv = 0.66;

areaB = mean(cvData(2).cv > threshCv);
areaC = mean(cvData(1).cv < threshCv);
areaA = mean(cvData(2).cv < threshCv) - areaC;

% NOTE(amotta): areas A + B + C = 1
learnLowerBound = areaA / (areaA + areaB + areaC) %#ok
learnUpperBound = (1 - areaB) / (areaA + areaB + areaC) %#ok

%% plot same-axon same-dendrite CV vs. intersynapse distance
maxInterSynDistUm = 25;

% calculate inter-synapse distance
dupNeurites.synAreaCv = ...
    std(dupNeurites.synAreas, 0, 2) ...
    ./ mean(dupNeurites.synAreas, 2);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

scatter(ax, ...
    dupNeurites.interSynDist / 1E3, ...
    dupNeurites.synAreaCv, 'x');

mask = (dupNeurites.interSynDist < maxInterSynDistUm * 1E3);
b = [ ...
    ones(size(dupNeurites.interSynDist(mask))), ...
    dupNeurites.interSynDist(mask) / 1E3] \ dupNeurites.synAreaCv(mask);

fitF = @(x) b(1) + x .* b(2);
rawName = sprintf('Raw data (n = %d)', numel(xLog));
fitName = sprintf( ...
    'Linear fit (y = %.2g + %.4gd; d_{max} = %g µm)', ...
    b(1), b(2), maxInterSynDistUm);

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));

ax.TickDir = 'out';
xlabel(ax, 'Intersynapse distance (µm)');
ylabel(ax, 'Coefficient of variation');
legend(ax, rawName, fitName, 'Location', 'SouthEast');

title( ...
   {'Distance dependence of synaptic consistency'; ...
    info.git_repos{1}.hash}, 'FontWeight', 'normal', 'FontSize', 10);

%% debugging
%{
%{
% select highly inconsistent pairs
randDupNeurites = find( ...
    dupNeurites.synAreas(:, 2) ...
 ./ dupNeurites.synAreas(:, 1) > 10);
%}

% select random pairs
randDupNeurites = 1:size(dupNeurites, 1);

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
%}