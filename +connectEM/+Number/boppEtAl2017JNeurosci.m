% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

%% TODO
% [ ] Take into account relative frequencies of VGluT2+/-

%% Configuration
% NOTE(amotta): This values were extracted from panel 5b of
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

numAxons = 100;
numSynsPerAxon = 10;

%% Simulate
rng(0);
ccVec = arrayfun(@(~) mean(simulateAxon( ...
    synPerBoutonProb(1, :), numSynsPerAxon) > 1), 1:numAxons);
tcVec = arrayfun(@(~) mean(simulateAxon( ...
    synPerBoutonProb(2, :), numSynsPerAxon) > 1), 1:numAxons);

%% Plot
figure; hold on;
binEdges = linspace(0, 1, 11);
histogram(ccVec, binEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(tcVec, binEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);

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
