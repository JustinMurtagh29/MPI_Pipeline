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

typeT = table;
typeT.id = { ...
    'Corticocortical'; 'Thalamocortical'; 'Inhibitory'; ... % axon classes
    'InhibitoryApicalDendrite'}; % specificity-based axon classes

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.synType = conn.axonMeta.axonClass(synT.preAggloId);

%% Specificity-based connectome
specConn = conn;
specAxonClasses = axonClasses(1:2);

[specConn, specAxonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        specConn, specAxonClasses, 'minSynPre', 10);
specAxonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses( ...
        specConn, specAxonClasses);

% Reset classes in axon meta data. The values are filled in below.
curSynTypes = {};
specConn.axonMeta.axonClass = ...
    repelem({[]}, numel(specConn.axons), 1);

% Modified from
% +connectEM/+Consistency/loadConnectome.m
% commit a23788a05313f61dfef152443cde5d7ff59d7622
for curIdx = 1:numel(specAxonClasses)
    curAxonClass = specAxonClasses(curIdx);
    curSpecs = curAxonClass.specs;
    
    curPrefix = strsplit(curAxonClass.title);
    curPrefix = curPrefix{1};
    curPrefix(1) = upper(curPrefix(1));
    
    curNames = fieldnames(curSpecs);
    curAxonIds = cellfun( ...
        @(n) curSpecs.(n).axonIds, ...
        curNames, 'UniformOutput', false);

    % Find axons with no and multiple specificities
   [curSpecAxonIds, ~, curDupAxonIds] = unique(cell2mat(curAxonIds));
    curDupAxonIds = curSpecAxonIds(accumarray(curDupAxonIds, 1) > 1);
    curNonSpecAxonIds = setdiff(curAxonClass.axonIds, curSpecAxonIds);

    curAxonIds = cellfun( ...
        @(ids) setdiff(ids, curDupAxonIds), ...
        curAxonIds, 'UniformOutput', false);
    
    curNames{end + 1} = 'MultiSpec'; %#ok
    curAxonIds{end + 1} = curDupAxonIds(:); %#ok
    
    curNames{end + 1} = 'NoSpec'; %#ok
    curAxonIds{end + 1} = curNonSpecAxonIds(:); %#ok
    
    for curNameIdx = 1:numel(curNames)
        curName = curNames{curNameIdx};
        curSpecs.(curName).axonIds = curAxonIds{curNameIdx};
    end
    
    % Update specificity classes
    specAxonClasses(curIdx).specs = curSpecs;
    
    % Update axon meta data
    curNames = strcat(curPrefix, curNames);
    curSynTypes = cat(1, curSynTypes(:), curNames(:));
    specConn.axonMeta.axonClass(cell2mat(curAxonIds)) = ...
        repelem(curNames, cellfun(@numel, curAxonIds), 1);
end

curSynTypes{end + 1} = 'Other';
specConn.axonMeta.axonClass(cellfun( ...
    @isempty, specConn.axonMeta.axonClass)) = {'Other'};

curSynTypes = categorical(curSynTypes, curSynTypes, 'Ordinal', true);
[~, curIds] = ismember(specConn.axonMeta.axonClass, curSynTypes);
specConn.axonMeta.axonClass = curSynTypes(curIds);

synT.specSynType = specConn.axonMeta.axonClass(synT.preAggloId);

%% Find synapses per type
clear cur*;
[~, curSynTypes] = ismember(synT.synType, typeT.id);
[~, curSpecSynTypes] = ismember(synT.specSynType, typeT.id);

typeT.synIds = accumarray( ...
   [nonzeros(curSynTypes(:)); nonzeros(curSpecSynTypes(:))], ...
   [find(curSynTypes(:)); find(curSpecSynTypes(:))], ...
   [height(typeT), 1], @(ids) {ids(:)}, {zeros(0, 1)});

%% Synapse positions
synPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = synPos(synT.id, :);
clear synPos;

% Physical units, relative to center
synT.pos = synT.pos - mean(param.bbox, 2)';
synT.pos = synT.pos .* param.raw.voxelSize ./ 1E3;

%% Numbers
clear cur*;
binUm = 5;

curRangeX = param.bbox(1, :) - mean(param.bbox(1, :), 2);
curRangeX = curRangeX .* param.raw.voxelSize(1) / 1E3;

% Correct for Benedikt's synapse margin
curRangeX = curRangeX + marginUm .* [+1, -1];

curBinEdges = [ ...
    -inf, curRangeX(1), curRangeX(1) + binUm, ...
    curRangeX(2) - binUm, curRangeX(2), +inf];

curSynData = nan(numel(curBinEdges) - 1, height(synT));
for curIdx = 1:height(typeT)
    curSynT = synT(typeT.synIds{curIdx}, :);
    curSynT.posId = curSynT.pos(:, 1);
    curSynT.posId = discretize(curSynT.posId, curBinEdges);
    
    curSynData(:, curIdx) = accumarray( ...
        curSynT.posId, 1, [numel(curBinEdges) - 1, 1]);
end

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

typeDensities = cell(height(typeT), 1);
for curIdx = 1:height(typeT)
    curSynT = synT(typeT.synIds{curIdx}, :);
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
typeBinCounts = nan(numel(typeBinEdges) - 1, height(typeT));

for curIdx = 1:height(typeT)
    curSynT = synT(typeT.synIds{curIdx}, :);
    curSynT.posId = curSynT.pos(:, dimIds(2));
    curSynT.posId = discretize(curSynT.posId, typeBinEdges);
    curSynT(isnan(curSynT.posId), :) = [];
    
    typeBinCounts(:, curIdx) = accumarray( ...
        curSynT.posId, 1, [numel(typeBinEdges) - 1, 1]);
end

%% Do the actual plotting
for curTypeIdx = 1:height(typeT)
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [600, 320];

    curIm = typeDensities{curTypeIdx};
    curIm = uint8(double(intmax('uint8')) * curIm / max(curIm(:)));
    curImSize = size(curIm);

    curImAx = axes(curFig); %#ok
    curIm = image(curImAx, curIm);
    colormap(curImAx, jet(256));
    axis(curImAx, 'image');

    curImAx.TickDir = 'out';
    
    curCortDim = find(dimIds == 1);
    curCortDim = char(double('X') - 1 + curCortDim);
    
    if isempty(curCortDim)
        curImAx.XTick = [];
        curImAx.YTick = [];
    else
        curImAx.([curCortDim, 'Tick']) = feval( ...
            @(v) v([1, end]), curImAx.([curCortDim, 'Tick']));
        curImAx.([setdiff('XY', curCortDim), 'Tick']) = [];
        curImAx.([curCortDim, 'TickLabel']) = {'Pia', 'WM'};
    end
    
    curSizeY = 0.6;
    curSizeX = curSizeY * curFig.Position(4) / curFig.Position(3);
    curSizeX = curSizeX * curImSize(2) / curImSize(1);
    curSize = [curSizeX, curSizeY];
    
    curImAx.Position(1:2) = (1 - curSize) / 2;
    curImAx.Position(3:4) = curSize;
    
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
    curHistAx.YAxis.Visible = 'off';

    curHistAx.Box = 'off';
    curHistAx.TickDir = 'out';
    curHistAx.YTick = [];
    curHistAx.YLim = typeBinEdges([1, end]);
    
    curPos = curImAx.Position;
    curHistAx.Position(1) = curPos(1) + curPos(3);
    curHistAx.Position(3) = 0.95 - curHistAx.Position(1);
    curHistAx.Position([2, 4]) = curPos([2, 4]);

    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            typeT.id{curTypeIdx}}, ...
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

curLeg = legend(curAx, typeT.id{curTypeIds}, 'Location', 'EastOutside');
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
