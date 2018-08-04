% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

% Ground truth annotations
oldConnFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% These axons were manually identified based on the following tracings:
% * TC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae5fff32700006943b3e2a8
% * CC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae6c9c9270000c253b3f60b
oldTcAxonIds = [ ...
    15197, 4391, 185, 6758, 5761, ...
    4793, 9960, 5712, 7195, 2162];
oldCcAxonIds = [ ...
    20115, 9052, 19696, 21957, 5040, ...
    13742, 7553, 28847, 20707, 18837];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

conn.axonMeta.fullPriSpineSynDens = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.pathLen;

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

boutonMetaFile = sprintf('%s_axonalBoutons_v1.mat', connName);
boutonMetaFile = fullfile(connDir, boutonMetaFile);
boutonMeta = load(boutonMetaFile);

%% Translate ground truth annotations to latest connectome
clear cur*;
curConn = load(oldConnFile);
curMaxSegId = Seg.Global.getMaxSegId(param);

curNewLUT = Agglo.buildLUT(curMaxSegId, conn.axons);
curFindInNew = @(segIds) mode(nonzeros(curNewLUT(segIds)));

tcAxonIds = cellfun(curFindInNew, curConn.axons(oldTcAxonIds));
tcAxonIds = reshape(tcAxonIds, size(oldTcAxonIds));

ccAxonIds = cellfun(curFindInNew, curConn.axons(oldCcAxonIds));
ccAxonIds = reshape(ccAxonIds, size(oldCcAxonIds));
clear cur*;

%% Build classifier
clear cur*;
curAxonIds = cat(1, ccAxonIds(:), tcAxonIds(:));
curY = cat(1, false(size(ccAxonIds(:))), true(size(tcAxonIds(:))));

curFeatNames = { ...
    'fullPriSpineSynDens', 'fullPriSpinesPerBouton', ...
    'fullPriSpinesMultiHitFrac', 'medianBoutonVol'};
curX = table2array(conn.axonMeta(:, curFeatNames));

classifier = fitclinear( ...
    curX(curAxonIds, :), curY, ...
    'Learner', 'logistic');

fprintf('\n');
fprintf('Results of logistic regression\n');
fprintf('* Bias: %f\n', classifier.Bias);
arrayfun(@(name, beta) fprintf( ...
    '* Weight for %s: %f\n', name{1}, beta), ...
    curFeatNames(:), classifier.Beta(:));

[~, axonProbs] = predict(classifier, curX);
axonProbs = axonProbs(:, 2);

%% Plot synapse densities for GT axons
clear cur*;
curConfigs = struct;
curConfigs(1).featName = 'fullPriSpineSynDens';
curConfigs(1).label = 'Spine synapse density (µm^{-1})';
curConfigs(1).binEdges = linspace(0, 0.4, 21);
curConfigs(1).thresh = 0.19;

curConfigs(2).featName = 'fullPriSpinesPerBouton';
curConfigs(2).label = 'Average number of spine synapses per bouton';
curConfigs(2).binEdges = linspace(0, 2.5, 26);
curConfigs(2).thresh = 1.55;

curConfigs(3).featName = 'fullPriSpinesMultiHitFrac';
curConfigs(3).label = 'Fraction of boutons with multiple spine synapses';
curConfigs(3).binEdges = linspace(0, 1, 21);
curConfigs(3).thresh = 0.475;

curConfigs(4).featName = 'medianBoutonVol';
curConfigs(4).label = 'Median bouton volume (µm^3)';
curConfigs(4).binEdges = linspace(0, 0.8, 21);
curConfigs(4).thresh = 0.3;

curExcAxonIds = axonClasses(1).axonIds;
curPlotProbs = true;

for curIdx = 1:numel(curConfigs)
    curConfig = curConfigs(curIdx);
    curBinEdges = curConfig.binEdges;
    curFeatVals = conn.axonMeta.(curConfig.featName);

    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [305, 285];

    curAx = axes(curFig); %#ok
    colors = curAx.ColorOrder;

    hold(curAx, 'on');

    yyaxis(curAx, 'right');
    curTcHist = histogram(curAx, curFeatVals(tcAxonIds));
    curCcHist = histogram(curAx, curFeatVals(ccAxonIds));
    
    ylim(curAx, [0, 1]);
    ylabel(curAx, 'Probability');

    yyaxis(curAx, 'left');
    curExcHist = histogram(curAx, curFeatVals(curExcAxonIds));
    curAllHists = [curTcHist, curCcHist, curExcHist];
    
    set(curAllHists, ...
        'DisplayStyle', 'stairs', ...
        'BinEdges', curBinEdges, ...
        'LineWidth', 2, 'FaceAlpha', 1);

   [curAllHists(1:2).Normalization] = deal('probability');
    curAllHists(1).EdgeColor = colors(1, :);
    curAllHists(2).EdgeColor = colors(2, :);
    curAllHists(3).EdgeColor = zeros(1, 3);
    
    if ~isempty(curConfig.thresh)
        line(curAx, ...
            repelem(curConfig.thresh, 2), curAx.YLim, ...
            'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');
    end

    curAx.TickDir = 'out';
    curAx.YColor = [0, 0, 0];
    xlim(curAx, curConfig.binEdges([1, end]));
    xlabel(curAx, curConfig.label);
    ylabel(curAx, 'Excitatory axons');

    curLeg = legend( ...
        [curTcHist, curCcHist], ...
        'Thalamocortical', ...
        'Corticocortical', ...
        'Location', 'SouthOutside');
    curLeg.Box = 'off';

    title( ...
        curAx, {info.filename, info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

%% Plot TC probability distribution
clear cur*;

% Configuration
curBinEdges = linspace(0, 1, 21);

curAxonIds = axonClasses(1).axonIds;
curTcProbs = conn.axonMeta.tcProb(curAxonIds);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [410, 195];

curAx = axes(curFig);
histogram(curAx, ...
    curTcProbs, ...
    'BinEdges', curBinEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2, ...
    'FaceAlpha', 1);

curAx.Box = 'off';
curAx.TickDir = 'out';
curAx.XLim = curBinEdges([1, end]);

xlabel(curAx, 'TC probability');
ylabel(curAx, 'Excitatory axons');
title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
