% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

minSynPre = 10;
info = Util.runInfo();

% Threshold p-values
% Chosen to be half of chance level.
pTheta = struct;
% Excitatory axons
pTheta(1).WholeCell = 0.02;
pTheta(1).ApicalDendrite = 0.035;
% Inhibitory axons
pTheta(2).Somata = 0.115;
pTheta(2).WholeCell = 0.1;
pTheta(2).ApicalDendrite = 0.06;
% Thalamocortical axons
pTheta(3).Somata = 0.015;
pTheta(3).WholeCell = 0.2;
pTheta(3).ApicalDendrite = 0.025;

%% Loading data
conn = ...
    connectEM.Connectome.load(param, connName);
[classConnectome, targetClasses] = ...
	connectEM.Connectome.buildClassConnectome(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

%% Perform analysis
axonClassIdx = 2;

% Generate list of target classes to check
targetClassIds = fieldnames(pTheta(axonClassIdx));
targetClassIds(cellfun(@(n) isempty(...
    pTheta(axonClassIdx).(n)), targetClassIds)) = [];
[~, targetClassIds] = ismember(targetClassIds, targetClasses);
targetClassIds = reshape(sort(targetClassIds), 1, []);

axonClass = axonClasses(axonClassIdx);
axonProbs = connectEM.Specificity.calcPoissonProbs( ...
    classConnectome, axonClass.axonIds, axonClass.poissAxonIds);

fig = figure();
for targetClassIdx = 1:numel(targetClassIds)
    targetClassId = targetClassIds(targetClassIdx);
    
    targetClassName = char(targetClasses(targetClassId));
    maxP = pTheta(axonClassIdx).(targetClassName);
    
    % Select specific axons
    specAxonIds = axonProbs(:, targetClassId);
    specAxonIds = axonClass.axonIds(specAxonIds < maxP);

    % Calculate conditional class connectome
    condTargetClassIds = setdiff( ...
        1:numel(targetClasses), targetClassId);
    condClassConnectome = ...
        connectEM.Connectome.buildClassConnectome( ...
            conn, targetClasses(condTargetClassIds));

    condClassConnectome = condClassConnectome(specAxonIds, :);
    condClassSynFrac = condClassConnectome ./ sum(condClassConnectome, 2);

    % Null hypothesis
    nullClassSynFrac = classConnectome( ...
        axonClass.poissAxonIds, condTargetClassIds);
    nullClassSynFrac = nullClassSynFrac ./ sum(nullClassSynFrac, 2);
    
    for curClassIdx = 1:numel(condTargetClassIds)
        curClassId = condTargetClassIds(curClassIdx);
        ax = subplot( ...
            numel(targetClassIds), numel(targetClasses), ...
            (targetClassIdx - 1) * numel(targetClasses) + curClassId);
        axis(ax, 'square');
        hold(ax, 'on');

        binEdges = linspace(0, 1, 11);
        histogram(ax, ...
            nullClassSynFrac(:, curClassIdx), ...
            'BinEdges', binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        histogram(ax, ...
            condClassSynFrac(:, curClassIdx), ...
            'BinEdges', binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        xlim(ax, [0, 1]);
        ylim(ax, [0, 1]);

        title( ...
            ax, char(targetClasses(curClassId)), ...
            'FontWeight', 'normal', 'FontSize', 10);
    end
    
    ax = subplot( ...
        numel(targetClassIds), numel(targetClasses), ...
        (targetClassIdx - 1) * numel(targetClasses) + targetClassId);
    ax.Visible = 'off';
    
    ann = annotation( ...
        'textbox', ax.Position, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'String', { ...
            sprintf('%s specific axon', targetClassName);
            sprintf('%d axons with p ≤ %g', numel(specAxonIds), maxP)});
end