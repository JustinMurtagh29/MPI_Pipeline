% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Perform analysis
axonIds = find(conn.axonMeta.synCount >= minSynPre);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

[synTypes, ~, synTypeIds] = unique(synT.type);

axonT = table;
axonT.id = reshape( ...
    1:numel(conn.axons), [], 1);
axonT.syns = accumarray( ...
   [synT.preAggloId, synTypeIds], ...
    1, [height(axonT), numel(synTypes)]);
axonT = axonT(axonIds, :);

%%
synCounts = sum(axonT.syns, 2);
synFracs = axonT.syns ./ synCounts;

%% preparations
nullSynTypeProbs = sum(axonT.syns, 1) / sum(axonT.syns(:));
axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
    axonT.syns, [], nullSynTypeProbs, 'distribution', 'binomial');

%% plotting
fig = figure;
fig.Color = 'white';
fig.Position(3:4) = [1850, 1025];

binEdges = linspace(0, 1, 21);
pValAxes = cell(size(synTypes));

for classIdx = 1:numel(synTypes)
    className = char(synTypes(classIdx));
    classProb = nullSynTypeProbs(classIdx);
    
    axonSynTypeFracs = synFracs(:, classIdx);
    axonClassNullProbs = axonNullProbs(:, classIdx);
    
    % Null hypothesis
    [nullSynFrac, nullAxonCount] = ...
        connectEM.Specificity.calcExpectedDist( ...
            synCounts, classProb, 'distribution', 'binomial');
    
    ksProb = ...
        connectEM.Specificity.kolmogorovSmirnovTest( ...
            axonSynTypeFracs, nullSynFrac, ...
            'nullWeights', nullAxonCount, ...
            'tail', 'smaller');
    
    nullBinId = discretize(nullSynFrac, binEdges);
    nullBinCount = accumarray(nullBinId, nullAxonCount);
    
    % Measured
    ax = subplot(3, numel(synTypes), classIdx);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        axonSynTypeFracs, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', nullBinCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
    xlabel(ax, 'Synapse fraction');
    ax.XAxis.TickDirection = 'out';
    ax.XAxis.Limits = [0, 1];
    
    ylabel(ax, 'Axons');
    ax.YAxis.TickDirection = 'out';
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YAxis.Scale = 'log';
        
    title(ax, ...
        {className; sprintf('p = %g (tailed KS)', ksProb)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    %% p-values
    curBinEdges = linspace(-1E-3, 1 + 1E-3, numel(binEdges));
    
   [expChanceProbs, expChanceCounts] = ...
        connectEM.Specificity.calcExpectedChanceProbDist( ...
            synCounts, classProb);
    
    curExpCounts = discretize(expChanceProbs, curBinEdges);
    curExpCounts = accumarray(curExpCounts, expChanceCounts);
    
    ax = subplot( ...
        3, numel(synTypes), ...
        numel(synTypes) + classIdx);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        axonClassNullProbs, ...
        'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    histogram(ax, ...
        'BinCounts', curExpCounts, ...
        'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
    ax.XLim = curBinEdges([1, end]);
    
    pValAxes{classIdx} = ax;
    
    %% alternative visualization
    % Compare p-value distribution against expectation:
    % We'd expect there to be `theta` percent of axons with a p-value
    % below `theta`. If there are, however, significantly more axons
    % with a p-value below `theta`, something interesting is going on.
    curPVal = sort(axonClassNullProbs, 'ascend');
    curPVal = reshape(curPVal, 1, []);
    
    ax = subplot( ...
        3, numel(synTypes), ...
        2 * numel(synTypes) + classIdx);
    axis(ax, 'square');
    hold(ax, 'on');
    
    [curPVal, ~, curPAxonFrac] = unique(curPVal);
    curPAxonFrac = accumarray(curPAxonFrac, 1);
    curPAxonFrac = cumsum(curPAxonFrac) / sum(curPAxonFrac);
    
    % Conservative estimate of false detection rate (FDR)
    curFdrEst = cumsum(expChanceCounts);
    curFdrEst = curFdrEst / curFdrEst(end);
    
    curFdrEst = interp1(expChanceProbs, curFdrEst, curPVal);
    curFdrEst = curFdrEst(:) ./ curPAxonFrac(:);
    
    curThetaIdx = 1 + find( ...
        curFdrEst(1:(end - 1)) <= 0.2 ...
        & curFdrEst(2:end) > 0.2, 1);
    
    plot(ax, curPVal, curFdrEst, 'LineWidth', 1);
    
    xlim(ax, [0, 1]);
    ylim(ax, [0, 1.2]);
    xlabel(ax, 'p-value');
    ylabel(ax, 'Estimated FDR');
    
    if ~isempty(curThetaIdx)
        curThetaPVal = curPVal(curThetaIdx);
        
        plot(ax, ...
            repelem(curThetaPVal, 2), ylim(ax), ...
            'Color', 'black', 'LineStyle', '--');
        title(ax, ...
            sprintf('p = %f', curThetaPVal), ...
            'FontWeight', 'normal', 'FontSize', 10);
    end
end
    
axes = fig.Children;
yMax = max(arrayfun(@(a) a.YLim(end), axes));
set(axes, 'YLim', [0, yMax]);
set(axes(1:3:end), 'YLim', [0, 1.2]);

% Legend
ax = axes(3);
axPos = ax.Position;
leg = legend(ax, ...
    'Observed', ...
    'Binomial model', ...
    'Location', 'East');
leg.Box = 'off';

% Fix positions
ax.Position = axPos;
leg.Position(1) = sum(axPos([1, 3])) + 0.005;

annotation( ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});
