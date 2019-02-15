% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%%
N = 5290;
cvKsPvalThresh = 2.11527e-31;
learnedFracs = [0.086, 0.379, 0.5];

mix = struct('mean', {}, 'std', {}, 'coeff', {});
mixNames = {};

% L4 → proximal dendrite synapses
% connectEM.Connectome.plotLayerLayerConnections 67e320ec88ac4369e8205118be2544f6dbc8de70
mix(1).mean = -0.70476; mix(1).std = 0.30209; mix(1).coeff = 1;
mixNames{1} = 'L4 → proximal dendrites synapses';

% L4 → apical dendrite synapses
% connectEM.Connectome.plotLayerLayerConnections 67e320ec88ac4369e8205118be2544f6dbc8de70
mix(2).mean = -0.90493; mix(2).std = 0.30629; mix(2).coeff = 1;
mixNames{2} = 'L4 → apical dendrites synapses';

% corticocortical primary spine synapses
% connectEM.Connectome.plotSynapseSizeConsistency 9476c84415afa274ce8e80df679014bb37172a72
mix(3).mean = -0.707029; mix(3).std = 0.298972; mix(3).coeff = 1;
mixNames{3} = 'Corticocortical primary spine synapses';

% NOTE(amotta): Mixing coefficient of the L4 connections
mixGrid = mix;
mixGrid(1).coeff = 10 .^ linspace(-3, 0, 61);
mixGrid(2).coeff = 10 .^ linspace(-3, 0, 61);

% Sanity check
assert(isequal(numel(mix), numel(mixNames)));
assert(isequal(fieldnames(mix), fieldnames(mixGrid)));
varNames = fieldnames(mix);

info = Util.runInfo();
Util.showRunInfo(info);

%% Grid search
clear cur*;
rng(0);

% Span grid
curGridVars = struct2cell(mixGrid);
grid = cell(1, 1, numel(curGridVars));
[grid{:}] = ndgrid(curGridVars{:});

grid = cat(1 + numel(curGridVars), grid{:});
grid = reshape(grid, [], numel(curGridVars));

%%
clear cur*;
rng(0);

fracs = nan(size(grid, 1), 3);
pVals = nan(size(grid, 1), 1);

% Precompute stuff
curRandn = randn(N, 2);
curRandperm = reshape(randperm(2 * N), N, 2);

for curId = 1:size(grid, 1)
    curMix = grid(curId, :);
    curMix = transpose(reshape(curMix, [], numel(mix)));
    curMix = array2table(curMix, 'VariableNames', varNames);
    
    curEdges = cumsum(curMix.coeff);
    curEdges = [0; curEdges(:) / curEdges(end)];
    
    curSaSdT = table;
    curSaSdT.gId = reshape(linspace(0, 1, N), [], 1);
    curSaSdT.gId = discretize(curSaSdT.gId, curEdges);

    curSaSdT.logAsi = curRandn .* curMix.std(curSaSdT.gId);
    curSaSdT.logAsi = curSaSdT.logAsi + curMix.mean(curSaSdT.gId);

    curSaSdT.asi = 10 .^ curSaSdT.logAsi;
    curSaSdT.cv = std(curSaSdT.asi, 0, 2) ./ mean(curSaSdT.asi, 2);

    curCtrlT = table;
    curCtrlT.asi = curSaSdT.asi(curRandperm);
    curCtrlT.cv = std(curCtrlT.asi, 0, 2) ./ mean(curCtrlT.asi, 2);
    
   [~, pVals(curId)] = kstest2( ...
        curSaSdT.cv, curCtrlT.cv, 'tail', 'larger');
   [curA, curB, curC] = ...
        calculateLearnedFraction(curSaSdT.cv, curCtrlT.cv);
    fracs(curId, :) = [curA, curB, curC];
    
    %% Plotting
    %{
    curBinEdges = linspace(0, 1.5, 31);
    
    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    histogram(curAx, ...
        curSaSdT.cv, 'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', 'LineWidth', 2);
    histogram(curAx, ...
        curCtrlT.cv, 'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', 'LineWidth', 2);
    
    curFig.Color = 'white';
    curFig.Position(3:4) = [280, 250];
    
    curAx.TickDir = 'out';
    xlim(curAx, curBinEdges([1, end]));
    xlabel(curAx, 'Coefficient of variation');
    ylabel(curAx, 'Occurences');
    
    curTitle = sprintf( ...
        'Parameter set #%d. %.1f %% learned', ...
        curId, 100 * fracs(curId, 1));
    title(curAx, { ...
        info.filename; info.git_repos{1}.hash; curTitle}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    %}
end

%% Plot "best" mix
clear cur*;
curGridVars = cellfun(@numel, struct2cell(mixGrid));
curGridVars = reshape(curGridVars, numel(varNames), numel(mix));

[curVarIds, curMixIds] = find(curGridVars > 1);
curDims = arrayfun(@(a, b) curGridVars(a, b), curVarIds, curMixIds);
assert(numel(curDims) == 2);

curPlots = struct;
curPlots(1).title = '-log_{10}(p-value)';
curPlots(1).data = -log10(pVals);
curPlots(2).title = 'Learned fraction';
curPlots(2).data = fracs(:, 1);

for curPlot = curPlots
    curIm = reshape(curPlot.data, curDims(1), curDims(2));

    curFig = figure();
    curAx = axes(curFig); %#ok
    imagesc(curAx, curIm);
    colorbar('peer', curAx);
    
    xlabel(curAx, sprintf( ...
        '%s of %s', varNames{curVarIds(2)}, mixNames{curMixIds(2)}));
    curAx.XTick = 1 + linspace(0, size(curIm, 2) - 1, 4);
    curAx.XTickLabel = arrayfun( ...
        @(id) num2str(feval(@(v) v(id), ...
            mixGrid(curMixIds(2)).(varNames{curVarIds(2)}))), ...
        curAx.XTick, 'UniformOutput', false);
    
    ylabel(curAx, sprintf( ...
        '%s of %s', varNames{curVarIds(1)}, mixNames{curMixIds(1)}));
    curAx.YTick = 1 + linspace(0, size(curIm, 1) - 1, 4);
    curAx.YTickLabel = arrayfun( ...
        @(id) num2str(feval(@(v) v(id), ...
            mixGrid(curMixIds(1)).(varNames{curVarIds(1)}))), ...
        curAx.YTick, 'UniformOutput', false);
    
    axis(curAx, 'square');
    curFig.Color = 'white';
    curAx.TickDir = 'out';
    curAx.YDir = 'normal';
    curAx.Box = 'off';
    
    title(curAx, ...
        {info.filename; info.git_repos{1}.hash; curPlot.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

%% Plot "best" mixture
clear cur*;

curMinId = abs(fracs(:, 1) - learnedFracs(1));
[~, curMinId] = min(curMinId);

curMix = grid(curMinId, :);
curMix = transpose(reshape(curMix, [], numel(mix)));
curMix = array2table(curMix, 'VariableNames', varNames);
curMix.weight = curMix.coeff / sum(curMix.coeff);

curLimits = [-2, 0.5];
curX = linspace(curLimits(1), curLimits(2), 101);
curY = nan(height(curMix), numel(curX));

curLegends = cell(1 + numel(mix), 1);
curLegends{end} = 'Mixture of Gaussians';

for curId = 1:height(curMix)
    curY(curId, :) = curMix.weight(curId) * ...
        normpdf(curX, curMix.mean(curId), curMix.std(curId));
    curLegends{curId} = sprintf( ...
        '%s (%.1f %% of mix)', ...
        mixNames{curId}, 100 * curMix.weight(curId));
end

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');
plot(curAx, curX, curY, 'LineWidth', 2);
plot(curAx, curX, sum(curY, 1), 'LineWidth', 2);

curFig.Color = 'white';
curAx.Box = 'off';
curAx.TickDir = 'out';

curLeg = legend(curAx, curLegends);
curLeg.Location = 'NorthWest';
curLeg.Box = 'off';

xlim(curAx, curX([1, end]));
xlabel(curAx, 'log_{10}(ASI area [µm²]');
ylabel(curAx, 'Probability density');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Utilities
function [learnedFrac, unlearnedFrac, cvThresh] = ...
        calculateLearnedFraction(pairCvs, ctrlCvs)
    t = table;
    t.cv = [pairCvs; ctrlCvs];
    t.pdf = [ ...
        +repelem(1 / numel(pairCvs), numel(pairCvs), 1); ...
        -repelem(1 / numel(ctrlCvs), numel(ctrlCvs), 1)];
    t = sortrows(t, 'cv', 'ascend');
    t.cdf = cumsum(t.pdf);
    
   [~, cvThreshIdx] = max(t.cdf);
    cvThresh = t.cv(cvThreshIdx);
    unlearnedFrac = mean(pairCvs > cvThresh);
    learnedFrac = t.cdf(cvThreshIdx);
end
