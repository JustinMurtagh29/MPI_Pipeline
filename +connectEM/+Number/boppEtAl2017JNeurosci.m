% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
rng(0);

%% Configuration
% NOTE(amotta): These following values were extracted from panel 5b of
% Bopp, Holler-Rickauer, Martin, Schuhknecht (2017). J Neurosci
%
% Rows corresponds to VGluT2- and VGluT2+
% Columns correspond to 1..4 synapses per axonal bouton
synPerBoutonProb = [ ...
    27.989,  4.828,  1.609, 0.000; ...
    11.638, 11.150, 10.172, 2.040];
synPerBoutonProb = ...
    synPerBoutonProb ...
 ./ sum(synPerBoutonProb, 2);

% NOTE(amotta): These following values were extracted from panel 6b of
% Bopp, Holler-Rickauer, Martin, Schuhknecht (2017). J Neurosci
%
% Rows corresponds to VGluT2- and VGluT2+
% Columns corresponds to bouton volumes.
boutonVolProb = transpose([ ...
    32.329,  4.770;
    30.978, 11.954;
     9.426,  5.977;
     8.075, 10.776;
     1.322,  3.592;
     0.000,  4.770;
     2.672,  7.184;
     0.000,  9.569;
     0.000,  8.362;
     1.322,  7.184;
     0.000,  4.770;
     0.000,  0.000;
     0.000,  2.385;
     0.000,  1.178;
     0.000,  0.000;
     0.000,  0.000;
     0.000,  1.178;
     0.000,  0.000;
     0.000,  1.178;
     0.000,  0.000;
     0.000,  1.178]);
boutonVolProb = ...
    boutonVolProb ...
 ./ sum(boutonVolProb, 2);

% See global "percentage of VGluT2 synapses" values for S1 in Table 1. This
% synapse fraction also roughly corresponds to the axon fraction under the
% following assumptions:
%
% * The following simulation is restricted to excitatory axons
% * Excitatory axons make exclusively asymmetric synapses
% * All excitatory axons make exactly the same number of synapses
tcSynFrac = 0.172;

numAxons = 5000;
numSynsPerAxon = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Simulate
tcMask = linspace(0, 1, numAxons) < tcSynFrac;

curSampleFun = @(probs) @(isTc) datasample( ...
    1:size(probs, 2), 1, 'Weights', probs(1 + isTc, :));

% Simulate individual boutons
boutonT = table;
boutonT.isTc = tcMask(:);
boutonT.numSyn = arrayfun(curSampleFun(synPerBoutonProb), boutonT.isTc);
boutonT.vol = arrayfun(curSampleFun(boutonVolProb), boutonT.isTc);

% Simulate entire axons
axonT = table;
axonT.isTc = tcMask(:);

axonT.synsPerBouton = arrayfun( ...
    @(isTc) simulateAxon( ...
        synPerBoutonProb(1 + isTc, :), numSynsPerAxon), ...
    axonT.isTc, 'UniformOutput', false);
axonT.boutonVols = arrayfun( ...
    @(isTc, numBoutons) datasample( ...
        1:size(boutonVolProb, 2), numBoutons, ...
        'Weights', boutonVolProb(1 + isTc, :)), ...
    axonT.isTc, cellfun(@numel, axonT.synsPerBouton), ...
    'UniformOutput', false);

axonT.avgNumSyn = cellfun(@mean, axonT.synsPerBouton);
axonT.fracMultiSyn = cellfun(@(n) mean(n > 1), axonT.synsPerBouton);
axonT.medianVol = cellfun(@median, axonT.boutonVols);

axonT(:, {'synsPerBouton', 'boutonVols'}) = [];

%% Plots for individual boutons
clear cur*;

curConfigs = struct;
curConfigs(1).table = boutonT;
curConfigs(1).prob = synPerBoutonProb;
curConfigs(1).var = 'numSyn';

curConfigs(2) = curConfigs(1);
curConfigs(2).prob = boutonVolProb;
curConfigs(2).var = 'vol';

curColors = get(groot, 'defaultAxesColorOrder');

for curConfig = curConfigs
    curIsTc = curConfig.table.isTc;
    curVar = curConfig.table.(curConfig.var);
    
    % Histogram
    curFig = figure();
    curFig.Color = 'white';
    curAx = axes(curFig); %#ok
    curAx.TickDir = 'out';
    
    hold(curAx, 'on');
    curBinEdges = (0:size(curConfig.prob, 2)) + 1 / 2;

    curHist = @(data, color) histogram(curAx, ...
        'BinCounts', histcounts(data, curBinEdges) / numAxons, ...
        'BinEdges', curBinEdges, 'EdgeColor', 'none', ...
        'FaceColor', color, 'FaceAlpha', 1);
    
    curHist(curVar, curColors(1, :));
    curHist(curVar(~curIsTc), 'black');
end

%% Plot precision / recall for boutons
% Individual features and then combination of them
clear cur*;

curTables = {boutonT, axonT};
for curTableIdx = 1:numel(curTables)
    curT = curTables{curTableIdx};
    curT = curT(randperm(height(curT)), :);
    curT.isTrain(:) = linspace(0, 1, height(curT)) < 1 / 2;

    curVars = curT.Properties.VariableNames;
    curVars = setdiff(curVars, {'isTc', 'isTrain'});

    curT.data = table2array(curT(:, curVars));
    curT.data = zscore(curT.data, 0, 1);
    
    % Logistic regression
    curModel = fitglm( ...
        curT.data(curT.isTrain, :), ...
        curT.isTc(curT.isTrain), ...
        'Distribution', 'binomial');

    curT.logReg = predict(curModel, curT.data);
    curVars{end + 1} = 'logReg'; %#ok
    curT = curT(~curT.isTrain, :);

    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    for curVarIdx = 1:numel(curVars)
       [curPrec, curRec, curThresh] = ...
           precisionRecall(curT.(curVars{curVarIdx}), curT.isTc);
        plot(curRec, curPrec, 'LineWidth', 2);
    end
    
    curF1Vec = 1 ./ (((1 ./ curPrec) + (1 ./ curRec)) / 2);
   [~, curMaxF1Idx] = max(curF1Vec);
   
    plot(curAx, ...
        curRec(curMaxF1Idx), curPrec(curMaxF1Idx), 'o', ...
        'Color', 'black', 'LineWidth', 2, 'MarkerSize', 10);

    curLeg = legend(curAx, curVars);
    set(curLeg, 'Location', 'SouthWest', 'Box', 'off');
    
    % Cosmetics
    curFig.Color = 'white';
    curFig.Position(3:4) = [189, 202];
    curAx.TickDir = 'out';
    
    grid(curAx, 'on');
    axis(curAx, 'square');
    xticks(curAx, 0:0.1:1); curAx.XTickLabel(2:2:end) = {''};
    yticks(curAx, 0:0.1:1); curAx.YTickLabel(2:2:end) = {''};
    
    xlabel(curAx, 'Recall');
    ylabel(curAx, 'Precision');
    title(curAx,  ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    fprintf('Maximum F1 Score\n');
    fprintf('* Precision: %.2f %%\n', 100 * curPrec(curMaxF1Idx));
    fprintf('* Recall: %.2f %%\n', 100 * curRec(curMaxF1Idx));
    fprintf('* Threshold: %.2f %%\n', 100 * curThresh(curMaxF1Idx));
    fprintf('\n');
end

%% Utilities
function synsPerBouton = simulateAxon(synPerBoutonProb, numSyn)
    synsPerBouton = datasample( ...
        1:numel(synPerBoutonProb), numSyn, ...
        'Weights', synPerBoutonProb);
    synsPerBouton = repelem(1:numSyn, synsPerBouton);
    
    % Simulate border effects
    off = randi(1 + numel(synsPerBouton) - numSyn);
    synsPerBouton = synsPerBouton(off:(off + numSyn - 1));
    
   [~, ~, synsPerBouton] = unique(synsPerBouton);
    synsPerBouton = accumarray(synsPerBouton, 1);
end

function [precVec, recVec, threshVec] = precisionRecall(data, labels)
    data = reshape(data, 1, []);
   [threshVec, sortIds] = sort(data, 'descend');
    label = reshape(labels(sortIds), 1, []);

    posAll = sum(labels);
    posPred = 1:numel(data);
    truePosPred = cumsum(label);

    precVec = truePosPred ./ posPred;
    recVec = truePosPred / posAll;
end
