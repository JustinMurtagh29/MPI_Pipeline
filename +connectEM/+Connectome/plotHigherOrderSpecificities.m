% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

minSynPre = 10;

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

info = Util.runInfo();

%% Loading data
conn = ...
    connectEM.Connectome.load(param, connName);
[classConnectome, targetClasses] = ...
	connectEM.Connectome.buildClassConnectome(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

%% Perform analysis
for axonClassIdx = 1:numel(axonClasses)
    plotAxonClass( ...
        info, classConnectome, targetClasses, ...
        axonClasses(axonClassIdx), pTheta(axonClassIdx));
end

%% Core function
function plotAxonClass(info, classConn, targetClasses, axonClass, pTheta)
    axonProbs = connectEM.Specificity.calcPoissonProbs( ...
        classConn, axonClass.axonIds, axonClass.poissAxonIds);
    
    % Find classes for which a threshold was set
    specClassIds = fieldnames(pTheta);
    specClassIds(cellfun(@(n) isempty(pTheta.(n)), specClassIds)) = [];
    
    % Convert class names to indices
   [~, specClassIds] = ismember(specClassIds, targetClasses);
    specClassIds = reshape(sort(specClassIds), 1, []);

    fig = figure();
    for specClassIdx = 1:numel(specClassIds)
        specClassId = specClassIds(specClassIdx);
        
        % Look up p-value threshold
        specClassName = char(targetClasses(specClassId));
        specThresh = pTheta.(specClassName);

        % Select specific axons
        specAxonIds = axonProbs(:, specClassId);
        specAxonIds = axonClass.axonIds(specAxonIds < specThresh);

        % Calculate conditional class connectome
        condClassIds = setdiff( ...
            1:numel(targetClasses), specClassId);
        condClassConn = classConn(specAxonIds, condClassIds);
        condClassSynFrac = condClassConn ./ sum(condClassConn, 2);

        % Null hypothesis
        nullClassSynFrac = classConn( ...
            axonClass.poissAxonIds, condClassIds);
        nullClassSynFrac = nullClassSynFrac ./ sum(nullClassSynFrac, 2);

        for curClassIdx = 1:numel(condClassIds)
            curClassId = condClassIds(curClassIdx);
            
            ax = subplot( ...
                numel(specClassIds), numel(targetClasses), ...
                (specClassIdx - 1) * numel(targetClasses) + curClassId);
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
            numel(specClassIds), numel(targetClasses), ...
            (specClassIdx - 1) * numel(targetClasses) + specClassId);
        hold(ax, 'on');
        
        % Generate invisible fake plot for legend
        plot(ax, nan, nan, 'LineWidth', 2);
        plot(ax, nan, nan, 'LineWidth', 2);
        
        legend(ax, ...
            'Expected (Poisson)', 'Observed', ...
            'Location', 'south');
        ax.Visible = 'off';
        
        annotation( ...
            'textbox', ax.Position, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'String', { ...
                sprintf('%s specific', specClassName);
                sprintf('%d axons with p ≤ %g', ...
                    numel(specAxonIds), specThresh)});
    end
    
    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'String', { ...
            axonClass.title;
            info.filename;
            info.git_repos{1}.hash});
        
end