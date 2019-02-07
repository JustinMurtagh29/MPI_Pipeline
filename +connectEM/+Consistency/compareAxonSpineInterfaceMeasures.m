% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, connFile] = ...
    connectEM.Consistency.loadConnectome(param);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph = graph(graph.borderIdx ~= 0, :);

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn, 'addBorderIdsVar', true);

%% Calculate ASI positions
curBorderMeta = fullfile(rootDir, 'globalBorder.mat');
curBorderMeta = load(curBorderMeta, 'borderSize', 'borderCoM');
curBorderMeta = structfun(@double, curBorderMeta, 'UniformOutput', false);

curWeightedMean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);

asiT.pos = cellfun( ...
    @(ids) curWeightedMean( ...
        curBorderMeta.borderSize(ids), ...
        curBorderMeta.borderCoM(ids, :)), ...
	asiT.borderIds, 'UniformOutput', false);
asiT.pos = round(cell2mat(asiT.pos));
clear cur*;

%% Calculate areas
rng(0);
randIds = randperm(height(asiT));

areas = ...
    connectEM.Consistency.axonSpineInterfaceArea( ...
        param, asiT(randIds(1:20), :), conn.axons, shAgglos);

%% Generate plots
labels = { ...
    'Helmstaedter et al. 2013 Nature', ...
    'de Vivo et al. 2017 Science', ...
    'Isosurface', 'Alpha Shape'};

for curA = 1:(size(areas, 2) - 1)
    for curB = (curA + 1):size(areas, 2)
        figure;
        hold on;
        axis equal;
        curFit = fitlm(areas(:, curA), areas(:, curB));
        scatter(areas(:, curA), areas(:, curB), '.');
        plot(xlim(), curFit.predict(xlim()'));

        legend({ ...
            sprintf('Data points (n = %d)', size(areas, 1)); ...
            sprintf('y = %.2g + %.2gx', curFit.Coefficients.Estimate)}, ...
            'Location', 'NorthWest');

        set(gcf, 'Color', 'white');
        set(gca, 'TickDir', 'out');

        xlim([0, 1.5]);
        xlabel({'ASI area à la'; labels{curA}});

        ylim([0, 1.5]);
        ylabel({'ASI area à la'; labels{curB}});
    end
end
