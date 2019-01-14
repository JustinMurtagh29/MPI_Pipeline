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

%% Plot
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

curFig.Position(3:4) = [390, 265];
curFig.Color = 'white';
curAx.TickDir = 'out';

xlim(curAx, curBinEdges([1, end]));
xlabel(curAx, 'Fraction of multi-synaptic boutons');
ylabel(curAx, 'Fraction of axons');
title(curAx,  ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

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
