% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

targetClasses = { ...
    'Somata', 'WholeCell', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;
info = Util.runInfo();

%% loading data
conn = load(connFile);

%% calculate target-class probability over synapses
[~, targetClassSyns] = ismember(conn.denMeta.targetClass, targetClasses);
targetClassSyns = accumarray(targetClassSyns, conn.denMeta.synCount);
targetClassProbs = targetClassSyns ./ sum(targetClassSyns);

%%
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

%% restrict to axons with enough synapses
axonIds = find(conn.axonMeta.synCount >= minSynPre);

[synCounts, ~, synCountAxons] = unique(conn.axonMeta.synCount(axonIds));
synCountAxons = accumarray(synCountAxons, 1);

specificities = classConnectome(axonIds, :);
specificities = specificities ./ sum(specificities, 2);

%% plotting
fig = figure;
binEdges = linspace(0, 1, 51);
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
    
    ax = subplot(1, numel(targetClasses), classIdx);
    hold(ax, 'on');
    
    histogram(ax, 'BinEdges', binEdges, 'BinCounts', poissBinCount, 'EdgeColor', 'none');
    histogram(ax, specificities(:, classIdx), binEdges, 'EdgeColor', 'none');

    xlabel(ax, sprintf('S(%s)', className));
    ax.XAxis.TickDirection = 'out';
    ax.XAxis.Limits = [0, 1];
    
    ax.YAxis.TickDirection = 'out';
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YAxis.Scale = 'log';
    
    if classIdx == 1
        legend('Expected (Poisson)', 'Observed');
    end
    
    axes{classIdx} = ax;
end

axes = horzcat(axes{:});
yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
for ax = axes; ax.YAxis.Limits(end) = yMax; end