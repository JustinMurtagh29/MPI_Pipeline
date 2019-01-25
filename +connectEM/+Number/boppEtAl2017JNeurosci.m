% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

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
%
% IMPORTANT: Note the transposition at the end!
boutonVolProb = [ ...
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
     0.000,  1.178]';
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

numAxons = 1000;
numSynsPerAxon = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Simulate
rng(0);

tcMask = linspace(0, 1, numAxons) < tcSynFrac;
fracMultiSynBoutons = arrayfun(@(isTc) mean(simulateAxon( ...
    synPerBoutonProb(1 + isTc, :), numSynsPerAxon) > 1), tcMask);
numSynsPerBouton = arrayfun(@(isTc) datasample( ...
    1:size(synPerBoutonProb, 2), 1, ...
    'Weights', synPerBoutonProb(1 + isTc, :)), tcMask);

%% Plot histogram
curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

curColors = get(groot, 'defaultAxesColorOrder');
curBinEdges = linspace(0, 1, 11);

curHist = @(data, color) histogram(curAx, ...
    'BinCounts', histcounts(data, curBinEdges) / numAxons, ...
    'BinEdges', curBinEdges, 'EdgeColor', 'none', ...
    'FaceColor', color, 'FaceAlpha', 1);
curHist(fracMultiSynBoutons, curColors(1, :));
curHist(fracMultiSynBoutons(~tcMask), 'black');

curLeg = legend('Thalamocortical axons', 'Corticocortical axons');
set(curLeg, 'Location', 'NorthEast', 'Box', 'off');

curFig.Position(3:4) = [261, 187];
curFig.Color = 'white';
curAx.TickDir = 'out';

xlim(curAx, curBinEdges([1, end]));
xlabel(curAx, 'Fraction of multi-synaptic boutons');
ylabel(curAx, 'Fraction of axons');
title(curAx,  ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot precision / recall
clear cur*;

[precVec, recVec, threshVec] = ...
    precisionRecall(fracMultiSynBoutons, tcMask);
[boutonPrecVec, boutonRecVec, boutonThreshVec] = ...
    precisionRecall(numSynsPerBouton, tcMask);

curF1Vec = 1 ./ (((1 ./ precVec) + (1 ./ recVec)) / 2);
[~, curMaxF1Idx] = max(curF1Vec);

curFig = figure();
curAx = axes(curFig);

hold(curAx, 'on');
grid(curAx, 'on');
axis(curAx, 'square');

curColors = get(groot, 'defaultAxesColorOrder');
plot(curAx, boutonRecVec, boutonPrecVec, 'LineWidth', 2, 'Color', 'black');
plot(curAx, recVec, precVec, 'LineWidth', 2, 'Color', curColors(1, :));

plot(curAx, ...
    recVec(curMaxF1Idx), precVec(curMaxF1Idx), 'o', ...
    'Color', 'black', 'LineWidth', 2, 'MarkerSize', 10);

curLeg = legend( ...
    'Individual axonal boutons', sprintf( ...
    'Axons with %d synapses', numSynsPerAxon));
set(curLeg, 'Location', 'SouthWest', 'Box', 'off');

curFig.Position(3:4) = [189, 202];
curFig.Color = 'white';
curAx.TickDir = 'out';

xticks(curAx, 0:0.1:1); curAx.XTickLabel(2:2:end) = {''};
yticks(curAx, 0:0.1:1); curAx.YTickLabel(2:2:end) = {''};

xlabel(curAx, 'Recall');
ylabel(curAx, 'Precision');
title(curAx,  ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

% Report numbers
fprintf('Maximum F1 Score\n');
fprintf('* Precision: %.2f %%\n', 100 * precVec(curMaxF1Idx));
fprintf('* Recall: %.2f %%\n', 100 * recVec(curMaxF1Idx));
fprintf('* Threshold: %.2f %%\n', 100 * threshVec(curMaxF1Idx));

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
   [threshVec, sortIds] = sort(data, 'descend');
    label = labels(sortIds);

    posAll = sum(labels);
    posPred = 1:numel(data);
    truePosPred = cumsum(label);

    precVec = truePosPred ./ posPred;
    recVec = truePosPred / posAll;
end
