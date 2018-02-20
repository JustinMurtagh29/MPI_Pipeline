% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

targetClasses = { ...
    'Somata', 'WholeCell', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;
maxSpineSynFracPoiss = 0.5;
info = Util.runInfo();

%% loading data
conn = load(connFile);

%% find inhibitory axons
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonIds = find(conn.axonMeta.synCount >= minSynPre);
poissAxonIds = find(conn.axonMeta.spineSynFrac < maxSpineSynFracPoiss);
poissAxonIds = intersect(axonIds, poissAxonIds);

%% build class connectome
% try to replicate cass connectome
connectome = conn.connectome;

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

[~, connectome.targetClassId] = ismember( ...
    connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);

%% calculate probabilities for Poisson
targetClassSyns = sum(classConnectome(poissAxonIds, :), 1);
targetClassProbs = targetClassSyns ./ sum(targetClassSyns);

%% restrict to axons with enough synapses
[synCounts, ~, synCountAxons] = unique(conn.axonMeta.synCount(axonIds));
synCountAxons = accumarray(synCountAxons, 1);

specificities = classConnectome ./ sum(classConnectome, 2);
specificities = specificities(axonIds, :);

%% plotting
fig = figure;
binEdges = linspace(0, 1, 26);
axes = cell(size(targetClasses));

for classIdx = 1:numel(targetClasses)
    className = targetClasses{classIdx};
    classProb = targetClassProbs(classIdx);
    
    % Poissons
    poiss = table;
    poiss.prob = cell2mat(arrayfun(@(nSyn, nAxons) ...
        nAxons * poisspdf((0:nSyn)', nSyn * classProb), ...
        synCounts, synCountAxons, 'UniformOutput', false));
    poiss.spec = cell2mat(arrayfun( ...
        @(nSyn) (0:nSyn)' ./ nSyn, ...
        synCounts, 'UniformOutput', false));
    
    poiss.binId = discretize(poiss.spec, binEdges);
    poissBinCount = accumarray(poiss.binId, poiss.prob);
    
    % Measured
    
    ax = subplot(1, numel(targetClasses), classIdx);
    hold(ax, 'on');
    
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', poissBinCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    histogram(ax, ...
        specificities(:, classIdx), binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);

    xlabel(ax, sprintf('S(%s)', className));
    ax.XAxis.TickDirection = 'out';
    ax.XAxis.Limits = [0, 1];
    
    ax.YAxis.TickDirection = 'out';
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YAxis.Scale = 'log';
    
    if classIdx == numel(targetClasses)
        legend(ax, ...
            'Expected (Poisson)', 'Observed', ...
            'Location', 'NorthEast');
    end
    
    axes{classIdx} = ax;
end

axes = horzcat(axes{:});
yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
for ax = axes; ax.YAxis.Limits(end) = yMax; end

annotation(...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', { ...
        'Observed specificities vs. Poisson model';
        sprintf('Based on axons with spine synapse fraction < %.1f', maxSpineSynFracPoiss);
        info.git_repos{1}.hash});