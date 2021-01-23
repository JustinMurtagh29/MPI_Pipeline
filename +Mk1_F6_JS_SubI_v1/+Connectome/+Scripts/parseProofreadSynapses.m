% After extract synaptic locations from connectome and connectomeMeta and generate nmls for manual proofreading
% use proofread nml to extract correct pairs and their post target types
%{
%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

[~, curConnName] = fileparts(connFile);
curVer = 'v1';

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
p = param.p;

conn = load(connFile);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'maxSegId');

% doing for sasd pairs
idxOut = arrayfun(@(x) numel(x{1})==2, conn.connectome.synIdx);
synIndices = conn.connectome.synIdx(idxOut); % maybe more than one synIdx per contact
synIndicesLookup = transpose(horzcat(synIndices{:}));
contactArea = conn.connectomeMeta.contactArea(idxOut);
contactAreaLookup = transpose(horzcat(contactArea{:}));
%}
%{
%% Specify which synapses file to parse after proofread
%nmlFile = fullfile(rootDir,'connectome','nmls','proofread',[curConnName,sprintf('-proofreadSynapses-%s-finished.nml', curVer)]);
% merged with proofread annotations
nmlFile = fullfile(rootDir,'connectome','nmls','proofread','MkL4-all-SASD-merged-finished.nml');
sprintf('Extracting sasd pairs in %s', nmlFile)

skel = skeleton(nmlFile);
treeNames = skel.names;

%% extract data per sasd pair
countId = 700;
sprintf('Searching pairs till Id %d', countId)
outTable = table;
outTable.idSASD = reshape(1:countId,'',1);
outTable.keep = nan(countId,1);
syn1 = {}; % type
syn2 = {}; % type
synIdx1 = [];
synIdx2 = [];
asi1 = [];
asi2 = [];

fun = @(x) str2double(x.id);

for idSASD = 1:countId
    curId = sprintf('sasd-%04d',idSASD);
    curPair = find(contains(treeNames, curId));
    assert(numel(curPair)==2)
    %if isempty(curPair)
    %    continue;
    %end
    curTree1 = treeNames{curPair(1)};
    curTree2 = treeNames{curPair(2)};

    idPre(idSASD) = fun(regexp(curTree1,'pre-(?<id>\d+)-','names'));
    idPost(idSASD) = fun(regexp(curTree1,'post-(?<id>\d+)-','names'));
    curSynIdx1 = fun(regexp(curTree1,'synIdx-(?<id>\d+)','names')); 
    curSynIdx2 = fun(regexp(curTree2,'synIdx-(?<id>\d+)','names'));

    % get contactArea from synIdx, Can be multiple matches so match row
    idxMatch = ismember(synIndicesLookup, [curSynIdx1, curSynIdx2], 'rows');
    flipFlag = false;

    if ~any(idxMatch)
        idxMatch = ismember(synIndicesLookup, [curSynIdx2, curSynIdx1], 'rows');
        flipFlag = true; % flipped syn1 and syn2
    end

    if ~any(idxMatch)
        error(sprintf('Could not find matching indices for the syns for idSASD %d', idSASD))
    end

    if flipFlag
        asiSynIdx1 = contactAreaLookup(idxMatch,2);
        asiSynIdx2 = contactAreaLookup(idxMatch,1);
    else
        asiSynIdx1 = contactAreaLookup(idxMatch,1);
        asiSynIdx2 = contactAreaLookup(idxMatch,2);
    end

    % save
    synIdx1(idSASD) = curSynIdx1;
    synIdx2(idSASD) = curSynIdx2;
    asi1(idSASD) = asiSynIdx1(1); % two pairs were duplicate 
    asi2(idSASD) = asiSynIdx2(1); 

    if contains(curTree1,'True') & ~contains(curTree1,{'merger','ign'})  & contains(curTree2,'True') & ~contains(curTree2,{'merger','ign'}) 
        outTable.keep(idSASD) = 1; % same-axon same-dendrite
        syn1{idSASD} = funType(regexp(curTree1,'(True-(?<type>\w+)','names'));
        syn2{idSASD} = funType(regexp(curTree2,'(True-(?<type>\w+)','names'));
    elseif (contains(curTree1,'True') & contains(curTree2,'True')) & (contains(curTree1,{'merger'}) | contains(curTree2,{'merger'})) % correct synapses but merger pre or post side
        outTable.keep(idSASD) = -1; % correct synapses but merger pre or post side
        syn1{idSASD} = funType(regexp(curTree1,'(True-(?<type>\w+)','names'));
        syn2{idSASD} = funType(regexp(curTree2,'(True-(?<type>\w+)','names'));
    else
        outTable.keep(idSASD) = 0; % false syn
        syn1{idSASD} = 'ign';
        syn2{idSASD} = 'ign';
    end
end

outTable.syn1 = reshape(syn1,'',1);
outTable.syn2 = reshape(syn2,'',1);
outTable.asi1 = reshape(asi1,'',1);
outTable.asi2 = reshape(asi2,'',1);

% keep only correct pairs in dataTable
dataTable = outTable(outTable.keep == 1,:); % keep true SASD

% keep pairs with mergers pre or post side in otherTable
otherTable = outTable(outTable.keep == -1,:); % keep non-SASD

%% Spine-Spine analysis: Asi-1 vs Asi2 Bartol et al 2015
% SASD pairs
idxPlot = contains(dataTable.syn1,{'Spine','Prim','Second'}) & contains(dataTable.syn2,{'Spine','Prim','Second'});
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;

d = sort([x1,x2],2,'descend'); % sort asi data
x1 = d(:,1); x2 = d(:,2);
clear d idxPlot
assert(numel(x1) == numel(x2))
sprintf('Found %d SASD spine-spine pairs',numel(x1))

% non-SASD pairs
idxPlot = contains(otherTable.syn1,{'Spine','Prim','Second'}) & contains(otherTable.syn2,{'Spine','Prim','Second'});
z1 = dataTable(idxPlot,:).asi1;
z2 = dataTable(idxPlot,:).asi2;

d = sort([z1,z2],2,'descend'); % sort asi data
z1 = d(:,1); z2 = d(:,2);
clear d idxPlot
assert(numel(z1) == numel(z2))
sprintf('Found %d non-SASD spine-spine pairs',numel(z1))
%}

%{
%% plot
fig = figure;
fig.Color = 'white';
ax = gca;
hold on
scatter(x1,x2,'kx');
grid('off');

curLimX = [1E-2, 1E1];
curLimY = [1E-2, 1E1];

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

ax.XAxis.Limits = curLimX;
ax.YAxis.Limits = curLimY;

xLog = log10(x1);
yLog = log10(x2);

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('y = %.2f x^{%.2f}', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));

plot(fitRange, fitRange, 'k--'); % diagonal

curLeg = legend(rawName, fitName);
set(curLeg, 'Box', 'Off', 'Location', 'best');

axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area [log_{10}(um^2)]')
ylabel('Asi 2 area [log_{10}(um^2)]')
set(gca,'FontSize',10)
title(ax, ...
    {info.filename; info.git_repos{1}.hash; 'Spine-Spine'}, ...
    'FontWeight', 'normal', 'FontSize', 10, 'Interpreter','none');
Util.setPlotDefault(ax)
outfile = fullfile(rootDir,'connectome','figures','scatter-Asi1-vs-Asi2-spine-spine-log-bartol_et_al.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all

%% Non log scale linear fit spine spine
curLimX = [0, +2];
curLimY = [0, +2];

fig = figure;
fig.Color = 'white';
ax = gca;
hold on

scatter(x1,x2,'kx');
grid('off')

curFit = fitlm(x1, x2);
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f) N=%d', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary, numel(x1)));
set(curLeg, 'Box', 'Off', 'Location', 'best');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area (um^2)')
ylabel('Asi 2 area (um^2)')
set(gca,'FontSize',10)
title(ax, ...
    {info.filename; info.git_repos{1}.hash; 'Spine-Spine'}, ...
    'FontWeight', 'normal', 'FontSize', 10, 'Interpreter','none');
Util.setPlotDefault(ax)
outfile = fullfile(rootDir,'connectome','figures','scatter-Asi1-vs-Asi2-spine-spine.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all
%}

%% histogram 1D
fig = figure;
fig.Color = 'white';
ax = gca;
hold on;
axis(ax, 'square');

cvData = cat(2, x1, x2);
curCvs = std(cvData, 0, 2) ./ mean(cvData, 2);

% Random pairs
rng(0);
randPairIds = randperm(2 * floor(numel(cvData)/ 2));
cvDataRandomized = cvData(reshape(randPairIds, [], 2));
% non SASD pairs
cvDataOther = cat(2, z1, z2);
% combined
%cvDataRand = cat(1, cvDataRandomized, cvDataOther);
cvDataRand  = cvDataRandomized;
curCvsRand = std(cvDataRand, 0, 2) ./ mean(cvDataRand, 2);

binEdges = linspace(0, 1.5, 16);

histogram( ...
     ax, curCvs, ...
    'BinEdges', binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

histogram( ...
     ax, curCvsRand, ...
    'BinEdges', binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

xlim(ax, binEdges([1, end]));
xlabel(ax, {'CV of synaspse size'});
ylabel(ax, 'Probability');

curLeg = legend(ax, sprintf('SASD (Mean ± std: %f ± %f)', mean(curCvs), std(curCvs)),...
                    sprintf('Random (Mean ± std: %f ± %f)', mean(curCvsRand), std(curCvsRand)));
set(curLeg, 'Box', 'Off', 'Location', 'Northeast');

title(ax, ...
    {info.filename; info.git_repos{1}.hash;}, ...
    'FontWeight', 'normal', 'FontSize', 10, 'Interpreter','none');
Util.setPlotDefault(ax)
outfile = fullfile(rootDir,'connectome','figures','variability_histogram_spine-spine.png')
export_fig(outfile,'-q101', '-nocrop','-m8')
close all


%% Mean synapse size vs Coeff. of variation
meanVec =  mean(cvData, 2);
cvVec = std(cvData, 0, 2) ./ mean(cvData, 2);

meanVecRand = mean(cvDataRand, 2);
cvVecRand = std(cvDataRand, 0, 2) ./ mean(cvDataRand, 2);

fig = figure; 
fig.Color = 'white';

subplot(1,2,1)
ax = gca;
hold on;
scatter(cvVec,log10(meanVec),'ko');

ax.XAxis.Limits = [0, 1.5];
ax.YAxis.Limits = [-1.5, 0.5];

ylabel('Mean synapse size [log_{10}(µm²)]');
xlabel({'CV of synapse sizes'});
title('Same-axon same-dendrite')
Util.setPlotDefault(ax)

subplot(1,2,2)
ax = gca;
hold on;
scatter(cvVecRand,log10(meanVecRand),'ko');

ax.XAxis.Limits = [0, 1.5];
ax.YAxis.Limits = [-1.5, 0.5];

ylabel('Mean synapse size [log_{10}(µm²)]');
xlabel({'CV of synapse sizes'});
title('Random pairs')
Util.setPlotDefault(ax)
%{
yLog = log10(meanVec); % do the fit
xLog = cvVec;

b = [ones(numel(yLog), 1), yLog] \ xLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('y = %.2f x^{%.2f}', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(yLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));

curLeg = legend(rawName, fitName);
set(curLeg, 'Box', 'Off', 'Location', 'best');
%}

h = suptitle(...
    {info.filename; info.git_repos{1}.hash});
set(h,'FontWeight', 'normal', 'FontSize', 10, 'Interpreter','none');

outfile = fullfile(rootDir,'connectome','figures','cv-vs-mean-asi-scatter-spine-spine.png')
export_fig(outfile,'-q101', '-nocrop','-m8')
close all

%% Heat map
%% Set p-value thresholds for evaluation
% NOTE(amotta): p-value threshold at which large and small low-CV regions
% are not yet merged. (I've only looked at range up to p = 10 %.)
curPvalThreshs = [0.5, 1, 2, 3, 4, 5] / 100;

% Excitatory axons
plotConfigs(1, 1).twoDimPValThreshs = 0.026;
%{
plotConfigs(2, 1).twoDimPValThreshs = inf;
plotConfigs(3, 1).twoDimPValThreshs = 0.015;
plotConfigs(4, 1).twoDimPValThreshs = inf;

% Corticocortical axons
plotConfigs(1, 2).twoDimPValThreshs = 0.052;
plotConfigs(2, 2).twoDimPValThreshs = inf;
plotConfigs(3, 2).twoDimPValThreshs = inf;
plotConfigs(4, 2).twoDimPValThreshs = inf;

% Thalamocortical axons
plotConfigs(1, 3).twoDimPValThreshs = 0.087;
plotConfigs(2, 3).twoDimPValThreshs = inf;
plotConfigs(3, 3).twoDimPValThreshs = inf;
plotConfigs(4, 3).twoDimPValThreshs = inf;
%}
for curIdx = 1:numel(plotConfigs)
    curThreshs = plotConfigs(curIdx).twoDimPValThreshs;
    
    % Only consider p-value thresholds before LTP and LTD regions merger
    curThreshs = curPvalThreshs(curPvalThreshs < curThreshs);
    
    % Only consider lowest and highest p-value threshold
    curThreshs = curThreshs(unique([1, numel(curThreshs)]));
    plotConfigs(curIdx).twoDimPValThreshs = curThreshs;
end


curConfig = struct;
curConfig.title = 'Same-axon same-dendrite';
curConfig.twoDimPValThreshs = plotConfigs(1,1).twoDimPValThreshs;

% Density difference map
curAxisX = 'cv';
curScaleY = 'log10';
curImSize = [256, 256];
curMethod = 'kde2d';

curRegionColors = [];
curRegionContourProps = {'LineColor', 'black'};

% NOTE(amotta): The `curMinMap` and `curMaxMap` matrices have the same size
% as the heatmaps of the CV × log10(avg. ASI) space. They contain the ASI
% areas of the smaller and larger synapses, respectively.
curLog10Avg = linspace(curLimY(1), curLimY(2), curImSize(1));
if strcmpi(curScaleY, 'linear'); curLog10Avg = log10(curLog10Avg); end
curCv = linspace(curLimX(1), curLimX(2), curImSize(2));
if strcmpi(curAxisX, 'reldiff'); curCv = curCv / sqrt(2); end
curCv(curCv >= sqrt(2)) = nan;

curMinMap = (10 .^ curLog10Avg(:)) .* (1 - curCv / sqrt(2));
curMaxMap = (10 .^ curLog10Avg(:)) .* (1 + curCv / sqrt(2));

% Show region in which the same-axon same-dendrite primary spine synapse
% pairs with at least one synapse in the smallest 10th percentile are
% located.
curMinPrctiles = 10;

switch lower(curAxisX)
    case 'cv'
        curLimX = [0, 1.5];
        curTicksX = linspace(0, 1.5, 4);
        curFuncX = @(a) std(a, 0, 2) ./ mean(a, 2);
        curLabelX = 'Coefficient of variation';
    case 'reldiff'
        curLimX = [0, 2];
        curTicksX = linspace(0, 2, 5);
        curFuncX = @(a) abs(diff(a, 1, 2)) ./ mean(a, 2);
        curLabelX = 'Relative difference';
    otherwise
        error('Invalid X axis "%s"', curAxisX);
end

switch lower(curScaleY)
    case 'linear'
        curLimY = [0, 1];
        curTicksY = linspace(0, 1, 5);
        curFuncY = @(areas) mean(areas, 2);
        curLabelY = 'Average ASI area [µm²]';
    case 'log10'
        curLimY = [-1.5, 0.5];
        curTicksY = linspace(-1.5, 0.5, 5);
        curFuncY = @(areas) log10(mean(areas, 2));
        curLabelY = 'log10(Average ASI area [µm²])';
    otherwise
        error('Invalid Y scale "%s"', curScaleY);
end

% SASD heatmap
curCtrlConfig = struct;
curCtrlConfig(1).title = 'SASD';

curKvPairs = { ...
     'xLim', curLimX, 'xAxis', curAxisX, ...
     'yLim', curLimY, 'yScale', curScaleY, ...
     'mapSize', curImSize, 'method', curMethod};
 
[curSaSdMap, curBw] = ...
        connectEM.Consistency.densityMap( ...
           cvData , curKvPairs{:});

curCtrlMaps = ...
    connectEM.Consistency.nullDensityMaps( ...
        cvData(:), curKvPairs{:}, ...
        'bandWidth', curBw, 'numMaps', 5000);

%% Prepare for figure
curCtrlMap = mean(curCtrlMaps, 3);

curMax = max(max(curSaSdMap(:)), max(curCtrlMap(:)));
curDiffMap = curSaSdMap - curCtrlMap;
curMaxDiff = max(abs(curDiffMap(:)));

curPvalMap = 1 - mean(curCtrlMaps < curSaSdMap, 3);
curPvalThreshs = sort(curConfig.twoDimPValThreshs, 'descend');

curPvalImg = -log10(min( ...
    1 - mean(curCtrlMaps < curSaSdMap, 3), ...
    1 - mean(curCtrlMaps > curSaSdMap, 3)));

%{
% NOTE(amotta): Detect statistically significant regions. Drop tiny
% regions (with less than 100 pixels or less than 1 % of SASD
% connections), which are most likely caused by outliers.
curRegionMask = curPvalMap < curPvalThreshs(1);
curRegionMask = curRegionMask & curRoiMask;
curRegionMask = bwlabel(curRegionMask);

curRegionProps = regionprops( ...
    curRegionMask, {'Area', 'Centroid', 'BoundingBox'}); %#ok
curKeepRegionIds = find([curRegionProps.Area] >= 100);

curKeepRegionIds = curKeepRegionIds(arrayfun( ...
    @(id) sum(curSaSdMap(curRegionMask(:) == id)), ...
    curKeepRegionIds) > 0.01);

curRegionProps = curRegionProps(curKeepRegionIds);
~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);
curSaSdT.regionId = curRegionMask(curSaSdT.mapIdx);

%}
curConfigTitle = sprintf( ...
    '%s (n = %d pairs)',...
    curConfig.title, size(cvData, 1));
curCtrlTitle = sprintf( ...
    'vs. random pairs of %s (n = %d)', ...
    curCtrlConfig.title,  floor(numel(cvData(:)) / 2) );

% Figure
curFig = figure();
curFig.Color = 'white';

curAx = subplot(2, 2, 1); % SASD map
imagesc(curAx, curSaSdMap);
caxis(curAx, [0, curMax]);
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
title('SASD map')

curAx = subplot(2, 2, 2); % Ctrl Map
imagesc(curAx, curCtrlMap);
caxis(curAx, [0, curMax]);
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
title('Randomized (Ctrl) map')

curAx = subplot(2, 2, 3); % p value map
curPValAx = curAx;

imagesc(curAx, curPvalImg);
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = arrayfun( ...
    @(val) sprintf('%.3g', val), ...
    curBar.Limits, 'UniformOutput', false);
curBar.Label.String = '-log10(p-value)';
title('-log10(p-value)')

%{
for curRegionId = 1:numel(curRegionProps)
    curPos = curRegionProps(curRegionId).Centroid;
    
    text(curAx, ...
        curPos(1), curPos(2), num2str(curRegionId), ...
        'Color', 'white', 'FontWeight', 'bold');
end
%}
curAx = subplot(2, 2, 4); % Diff map
hold(curAx, 'on');

imagesc(curAx, curDiffMap);
caxis(curAx, [-1, +1] * curMaxDiff);
colormap(curAx, jet(256));
%{
for curPvalThresh = curPvalThreshs
    contour(curAx, ...
        curRegionMask & ...
        curPvalMap < curPvalThresh, ...
        true, 'LineColor', 'black');
end
%}
curBar = colorbar('peer', curAx);
curBar.Ticks = [ ...
    curBar.Limits(1), ...
    mean(curBar.Limits), ...
    curBar.Limits(end)];
curBar.TickLabels = { ...
    sprintf('%.3g', -curMaxDiff), '0', ...
    sprintf('%.3g', +curMaxDiff)};
title('Diff map')
 
set( ...
       findobj(curFig.Children, 'Type', 'ColorBar'), ...
       'Location', 'EastOutside', ...
       'TickDirection', 'out', ...
       'Box', 'off');

curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
arrayfun(@(ax) hold(ax, 'on'), curAxes);

curTickIdsX = 1 + floor((curImSize(2) - 1) * ...
     (curTicksX - curLimX(1)) / (curLimX(2) - curLimX(1)));
curTickLabelsX = arrayfun( ...
    @num2str, curTicksX, 'UniformOutput', false);

curTickIdsY = 1 + floor((curImSize(1) - 1) * ...
    (curTicksY - curLimY(1)) / (curLimY(2) - curLimY(1)));
curTickLabelsY = arrayfun( ...
    @num2str, curTicksY, 'UniformOutput', false);

set(curAxes, ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'YDir', 'normal', ...
            'YTick', curTickIdsY, ...
            'YTickLabels', curTickLabelsY, ...
            'YLim', [1, curImSize(2)], ...
            'XTick', curTickIdsX, ...
            'XTickLabels', curTickLabelsX, ...
            'XLim', [1, curImSize(1)], ...
            'PlotBoxAspectRatio', [1, 1, 1], ...
            'DataAspectRatioMode', 'auto');      

arrayfun(@(ax) xlabel(ax, curLabelX), curAxes);
arrayfun(@(ax) ylabel(ax, curLabelY), curAxes);        

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', { ...%info.filename; info.git_repos{1}.hash; ...
        curConfigTitle; curCtrlTitle}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize',8);

outfile = fullfile(rootDir,'connectome','figures','sasd-heatmap_spine-spine.png')
export_fig(outfile,'-q101', '-nocrop','-m8')
close all

%{
% 2. plot shaft-shaft
fig = figure;
fig.Color = 'white';
ax = gca;
hold on;
idxPlot = contains(dataTable.syn1,{'Shaft'}) & contains(dataTable.syn2,{'Shaft'});
x1 = dataTable(idxPlot,:).asi1;
x2 = dataTable(idxPlot,:).asi2;
% sort data
d = sort([x1,x2],2,'descend');
x1 = d(:,1); x2 = d(:,2);
clear d

scatter(x1,x2,'kx')
count = sum(idxPlot);

curFit = fitlm(x1, x2);
plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f) N=%d', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary,count));
set(curLeg, 'Box', 'Off', 'Location', 'South');

ax.YAxis.Limits = curLimX;
ax.XAxis.Limits = curLimY;
axis('square')
ax.LineWidth = 2;
xlabel('Asi 1 area (um^2)')
ylabel('Asi 2 area (um^2)')
set(gca,'FontSize',10)
title(ax, ...
    {info.filename; info.git_repos{1}.hash; 'Shaft-Shaft'}, ...
    'FontWeight', 'normal', 'FontSize', 10);
outfile = fullfile(rootDir,'connectome','figures','sasd-pairs-ASI-shaft-shaft.png')
export_fig(outfile,'-q101', '-nocrop','-transparent')
close all

%% statistics
sprintf('False: %d', sum(~outTable.keep))
sprintf('Spine-Spine: %d', sum(contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Spine')))
sprintf('Spine-Prim: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Prim')) | ...
                                (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Spine')) ))
sprintf('Spine-Second: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Spine')) ))
sprintf('Prim-Prim: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Prim')) | ...
                                (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Prim')) ))

sprintf('Second-Second: %d', sum( (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Second')) ))

sprintf('Prim-Second: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Second')) | ...
                                (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Prim')) ))

sprintf('Spine-Shaft: %d', sum( (contains(dataTable.syn1,'Spine') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Spine')) ))
sprintf('Prim-Shaft: %d', sum( (contains(dataTable.syn1,'Prim') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Prim')) ))
sprintf('Second-Shaft: %d', sum( (contains(dataTable.syn1,'Second') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Second')) ))
sprintf('Shaft-Shaft: %d', sum( (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Shaft')) ))

sprintf('Soma-Soma: %d', sum( (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Soma')) | ...
                                (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Soma')) ))

sprintf('Soma-Shaft: %d', sum( (contains(dataTable.syn1,'Soma') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Soma')) ))

sprintf('Neck-Shaft: %d', sum( (contains(dataTable.syn1,'Neck') & contains(dataTable.syn2,'Shaft')) | ...
                                (contains(dataTable.syn1,'Shaft') & contains(dataTable.syn2,'Neck')) ))

%}
function type = funType(x)
    if isempty(x)
        type = 'Spine';
    else
        type = x.type;
    end
end

