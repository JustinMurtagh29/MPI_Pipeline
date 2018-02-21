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
info = Util.runInfo();

%% loading data
conn = load(connFile);

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

%% build axon classes
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonClasses = struct;
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.7);
axonClasses(1).title = 'axons with spine synapse fraction > 0.7';

axonClasses(2).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.5);
axonClasses(2).title = 'axons with spine synapse fraction > 0.5';

axonClasses(3).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac < 0.5);
axonClasses(3).title = 'axons with spine synapse fraction < 0.5';

axonClasses(4).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac < 0.3);
axonClasses(4).title = 'axons with spine synapse fraction < 0.3';

%% plot
for curIdx = 1:numel(axonClasses)
    plotAxonClass( ...
        info, conn.axonMeta, classConnectome, ...
        targetClasses, axonClasses(curIdx));
end

%% look at AD-specific axons
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

%%
className = 'ApicalDendrite';

axonSpecs = ...
    classConnectome ...
    ./ sum(classConnectome, 2);
axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & abs(conn.axonMeta.spineSynFrac - 0.5) > 0.2 ...
  & axonSpecs(:, strcmpi(targetClasses, className)) > 0.45);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', mfilename, info.git_repos{1}.hash));

% axons
for curIdx = 1:numel(axonIds)
    curAxonId = axonIds(curIdx);
    curAgglo = conn.axons{curAxonId};
    
    skel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    skel.names{end} = sprintf( ...
        'Axon %d. Spine fraction %.2f', ...
        curAxonId, conn.axonMeta.spineSynFrac(curAxonId));
    skel.colors{end} = [1, 0, 0, 1];
end

% targets
dendIds = find(conn.denMeta.targetClass == className);
for curIdx = 1:numel(dendIds)
    curDendId = dendIds(curIdx);
    curAgglo = conn.dendrites{curDendId};
    
    skel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    skel.names{end} = sprintf('Dendrite %d', curDendId);
    skel.colors{end} = [0, 0, 1, 1];
end

skel.write('/home/amotta/Desktop/test.nml');

%% plotting
function plotAxonClass(info, axonMeta, classConn, targetClasses, axonClass)
    %% preparations
    % calculate probabilities for Poisson
    targetClassSyns = sum(classConn(axonClass.axonIds, :), 1);
    targetClassProbs = targetClassSyns / sum(targetClassSyns);
    
    specificities = classConn(axonClass.axonIds, :);
    specificities = specificities ./ sum(specificities, 2);

    % restrict to axons with enough synapses
   [synCounts, ~, synCountAxons] = unique( ...
        axonMeta.synCount(axonClass.axonIds));
    synCountAxons = accumarray(synCountAxons, 1);
    
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
        'String', ...
           {'Observed specificities vs. Poisson model'; ...
            axonClass.title; info.git_repos{1}.hash});
end