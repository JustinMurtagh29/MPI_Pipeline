% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

dimIds = [3, 1];
maxImSize = 300;
bandWidth = [2, 2];
binSizeUm = 2;
marginUm = 3;

synTypes = { ...
    'Corticocortical', ...
    'Thalamocortical', ...
    'Inhibitory'};

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Preparing data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.synType = conn.axonMeta.axonClass(synT.preAggloId);

[~, synT.synType] = ismember(synT.synType, synTypes);
synT(~synT.synType, :) = [];

synPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = synPos(synT.id, :);
clear synPos;

% Physical units, relative to center
synT.pos = synT.pos - mean(param.bbox, 2)';
synT.pos = synT.pos .* param.raw.voxelSize ./ 1E3;

%% Numbers
binUm = 5;

curRangeX = param.bbox(1, :) - mean(param.bbox(1, :), 2);
curRangeX = curRangeX .* param.raw.voxelSize(1) / 1E3;

% Correct for Benedikt's synapse margin
curRangeX = curRangeX + marginUm .* [+1, -1];

curBinEdges = [ ...
    -inf, curRangeX(1), curRangeX(1) + binUm, ...
    curRangeX(2) - binUm, curRangeX(2), +inf];
curSynData = accumarray(horzcat( ...
    discretize(synT.pos(:, 1), curBinEdges), synT.synType), ...
    1, [numel(curBinEdges) - 1, numel(synTypes)]);

curSynData = curSynData([1 + 1, end - 1], :);

curInhCount = curSynData(:, 3);
curTcCount = curSynData(:, 2);

curInhExcRatio = curSynData(:, 3) ./ sum(curSynData(:, 1:3), 2);
curTcCcRatio = curSynData(:, 2) ./ sum(curSynData(:, 1:2), 2);

fprintf('Top %d µm of dataset\n', binUm);
fprintf('* Inh synapses: %d\n', curInhCount(1));
fprintf('* TC synapses: %d\n', curTcCount(1));
fprintf('* Inh / (Inh + Exc): %f\n', curInhExcRatio(1));
fprintf('* TC / (TC + CC): %f\n', curTcCcRatio(1));
fprintf('\n');

fprintf('Bottom %d µm of dataset\n', binUm);
fprintf('* Inh synapses: %d\n', curInhCount(end));
fprintf('* TC synapses: %d\n', curTcCount(end));
fprintf('* Inh / (Inh + Exc): %f\n', curInhExcRatio(end));
fprintf('* TC / (TC + CC): %f\n', curTcCcRatio(end));
fprintf('\n');

%% Scatter plot of synpapse position
clear cur*;
halfBoxUm = param.bbox(:, 2) - mean(param.bbox, 2);
halfBoxUm = halfBoxUm .* param.raw.voxelSize(:) / 1E3;
halfBoxUm = halfBoxUm - marginUm;

limits = binSizeUm * floor(halfBoxUm / binSizeUm);
limX = [-1, +1] .* limits(dimIds(1));
limY = [-1, +1] .* limits(dimIds(2));

%% Precompute images for faster prototyping
clear cur;
curImSize = [limY(2), limX(2)];
curImSize = round(maxImSize * curImSize / max(curImSize));

[curImGridY, curImGridX] = ndgrid( ...
    linspace(limY(1), limY(2), curImSize(1)), ...
    linspace(limX(1), limX(2), curImSize(2)));
curImGrid = cat(2, curImGridY(:), curImGridX(:));

typeDensities = cell(size(synTypes));
for curIdx = 1:numel(synTypes)
    curSynT = synT(synT.synType == curIdx, :);
    curSynT.posX = curSynT.pos(:, dimIds(1));
    curSynT.posY = curSynT.pos(:, dimIds(2));
    
    % Restrict to synapses in domain
    curSynT(abs(curSynT.posX) > limX(end), :) = [];
    curSynT(abs(curSynT.posY) > limY(end), :) = [];
    
    curIm = ksdensity( ...
       cat(2, curSynT.posY, curSynT.posX), curImGrid, ...
       'Support', cat(2, limY(:), limX(:)), ...
       'BoundaryCorrection', 'reflection', ...
       'BandWidth', bandWidth);
    curIm = reshape(curIm, curImSize);
    
    typeDensities{curIdx} = curIm;
end

%% Synapse histogram along Y axis of plot
clear cur*;
typeBinEdges = limY(1):binSizeUm:limY(2);

curSynT = synT;
curSynT.posId = curSynT.pos(:, dimIds(2));
curSynT.posId = discretize(curSynT.posId, typeBinEdges);
curSynT(isnan(curSynT.posId), :) = [];

typeBinCounts = accumarray( ...
    cat(2, curSynT.posId, curSynT.synType), ...
	1, [numel(typeBinEdges) - 1, numel(synTypes)]);

%% Do the actual plotting
for curTypeIdx = 1:numel(synTypes)
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [600, 320];
    curFigRatio = curFig.Position(4) / curFig.Position(3);

    curIm = typeDensities{curTypeIdx};
    curIm = uint8(double(intmax('uint8')) * curIm / max(curIm(:)));
    curImRatio = size(curIm, 1) / size(curIm, 2);

    curImAx = axes(curFig); %#ok
    curIm = image(curImAx, curIm);
    colormap(curImAx, jet(256));
    axis(curImAx, 'image');

    curImAx.XTick = [];
    curImAx.YTick = curImAx.YTick([1, end]);
    curImAx.YTickLabel = {'Pia', 'WM'};
    
    curImAx.Position(3:4) = 0.6 * [1, curImRatio / curFigRatio];
    curImAx.Position(1:2) = (1 - curImAx.Position(3:4)) / 2;
    
    % Linear regression
    curLinFit = (typeBinEdges(1:(end - 1)) + typeBinEdges(2:end)) / 2;
    curLinFit = fit(curLinFit(:), typeBinCounts(:, curTypeIdx), 'poly1');

    curHistAx = axes(curFig); %#ok
    hold(curHistAx, 'on');
    
    curHist = histogram(curHistAx, ...
        'BinEdges', typeBinEdges, ...
        'BinCounts', typeBinCounts(:, curTypeIdx), ...
        'DisplayStyle', 'stairs', 'LineWidth', 2, ...
        'FaceAlpha', 1, 'Orientation', 'horizontal');
    curHistFit = plot(curHistAx, ...
        curLinFit(typeBinEdges([1, end])), typeBinEdges([1, end]), ...
        'Color', 'black', 'LineWidth', 2);
    curHistAx.YDir = 'reverse';

    curHistAx.Box = 'off';
    curHistAx.TickDir = 'out';
    curHistAx.YTick = [];
    curHistAx.YLim = typeBinEdges([1, end]);
    
    curPos = curImAx.Position;
    curHistAx.Position([2, 4]) = curPos([2, 4]);
    curHistAx.Position(1) = curPos(1) + curPos(3);
    curHistAx.Position(3) = 0.95 - curHistAx.Position(1);
    curHistAx.Position(4) = curPos(4);

    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            synTypes{curTypeIdx}}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

%% Synapse density histogram
clear cur*;
curTypeIds = [1, 3];

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [380, 250];

curAx = axes(curFig);
axis(curAx, 'square');
hold(curAx, 'on');
for curTypeId = curTypeIds
    curColor = curAx.ColorOrder(curTypeId, :);
    
    histogram(curAx, ...
        'BinEdges', typeBinEdges, ...
        'BinCounts', typeBinCounts(:, curTypeId), ...
        'Orientation', 'horizontal', ...
        'DisplayStyle', 'stairs', ...
        'EdgeColor', curColor, ...
        'LineWidth', 2, ...
        'FaceAlpha', 1)
end

xlabel(curAx, 'Synapses per volume');
ylabel(curAx, 'Cortical depth');

curLimX = max(max(typeBinCounts(:, curTypeIds)));
curLimX = [0, 1E3 * ceil(curLimX / 1E3)];
curAx.XLim = curLimX;
curAx.YLim = limY;

curAx.YAxis.Direction = 'reverse';
curAx.TickDir = 'out';

curLeg = legend(curAx, synTypes(curTypeIds), 'Location', 'EastOutside');
curLeg.Box = 'off';

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Synapse ratios
clear cur*;

curConfigs = struct;
curConfigs(1).name = 'I / (I + E)';
curConfigs(1).binCounts = ...
    typeBinCounts(:, 3) ...
 ./ sum(typeBinCounts, 2);

curConfigs(2).name = 'TC / (TC + CC)';
curConfigs(2).binCounts = ...
    typeBinCounts(:, 2) ...
 ./ sum(typeBinCounts(:, 1:2), 2);

curX = (typeBinEdges(1:(end - 1)) + typeBinEdges(2:end)) / 2;
curLinFit = fit(curX(:), curConfigs(2).binCounts, 'poly1');

fitlm(curX(:), curConfigs(2).binCounts)
fprintf('* Number of CC synapses: %d\n', sum(typeBinCounts(:, 1)));
fprintf('* Number of TC synapses: %d\n', sum(typeBinCounts(:, 2)));
fprintf('* Number of synapses: %d\n', sum(sum(typeBinCounts(:, 1:2))));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [380, 250];

curAx = axes(curFig);
axis(curAx, 'square');
hold(curAx, 'on');

for curConfig = curConfigs
    histogram(curAx, ...
        'BinEdges', typeBinEdges, ...
        'BinCounts', curConfig.binCounts, ...
        'Orientation', 'horizontal', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1)
end

plot(curAx, ...
    curLinFit(typeBinEdges([1, end])), typeBinEdges([1, end]), ...
    'LineWidth', 2, 'Color', 'black');

xlabel(curAx, 'Synapse ratio per volume');
ylabel(curAx, 'Cortical depth');

curAx.YLim = limY;
curAx.YAxis.Direction = 'reverse';
curAx.TickDir = 'out';

curLeg = legend(curAx, {curConfigs.name}, 'Location', 'EastOutside');
curLeg.Box = 'off';

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
