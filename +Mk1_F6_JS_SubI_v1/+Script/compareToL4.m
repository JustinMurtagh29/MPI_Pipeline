% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
l4Config = struct;
l4Config.rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
l4Config.shFile = fullfile(l4Config.rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
l4Config.asiFile = fullfile(l4Config.rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified__20190227T082543_asiT.mat');

l23Config = struct;
l23Config.rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
l23Config.shFile = '/tmpscratch/amotta/l23/2019-11-18-axon-and-dendrite-agglomerates/20191203T021242_results_auto-spines_v2.mat';
l23Config.synTypeFile = fullfile(l23Config.rootDir, 'connectome', 'SynapseAgglomerates__20191203T021242_results__20191203T021242_results_auto-spines_v2__v1__types_v1.mat');
l23Config.asiFile = fullfile(l23Config.rootDir, 'connectome', 'Connectome_20191203T021242-results_20191203T021242-results-auto-spines-v2_SynapseAgglomerates--20191203T021242-results--20191203T021242-results-auto-spines-v2--v1__20191213T143516_asiT.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
% Layer 4
l4Config.param = Util.load(fullfile( ...
    l4Config.rootDir, 'allParameter.mat'), 'p');
l4ShAgglos = Util.load(l4Config.shFile, 'shAgglos');

l4SegSizes = Seg.Global.getSegToSizeMap(l4Config.param);
l4ShVols = cellfun(@(ids) sum(l4SegSizes(ids)), l4ShAgglos);
l4ShVols = l4ShVols * prod(l4Config.param.raw.voxelSize / 1E3);
clear l4ShAgglos l4SegSizes;

l4AsiT = Util.load(l4Config.asiFile, 'asiT');
l4AsiT = connectEM.Consistency.Calibration.apply(l4AsiT);
l4AsiT = l4AsiT(l4AsiT.type == 'PrimarySpine', :);
l4AsiT = l4AsiT(l4AsiT.area > 0, :);

% Layer 2/3
l23Config.param = Util.load(fullfile( ...
    l23Config.rootDir, 'allParameter.mat'), 'p');
[l23ShAgglos, l23ShAttached] = Util.load( ...
    l23Config.shFile, 'shAgglos', 'attached');

l23SegSizes = Seg.Global.getSegToSizeMap(l23Config.param);
l23ShVols = cellfun(@(ids) sum(l23SegSizes(ids)), l23ShAgglos);
l23ShVols = l23ShVols * prod(l23Config.param.raw.voxelSize / 1E3);
clear l23ShAgglos l23SegSizes;

l23AsiT = Util.load(l23Config.asiFile, 'asiT');
curTypes = Util.load(l23Config.synTypeFile, 'types');
l23AsiT.type = curTypes(l23AsiT.id);
l23AsiT = connectEM.Consistency.Calibration.apply(l23AsiT);
l23AsiT = l23AsiT(l23AsiT.type == 'PrimarySpine', :);
l23AsiT = l23AsiT(l23AsiT.area > 0, :);

%% Plot spine head volumes
clear cur*;
curBinEdges = linspace(-3, 0, 41);

fig = figure();
ax = axes(fig);

hold(ax, 'on');
histogram(ax, log10(l4ShVols), 'BinEdges', curBinEdges);
histogram(ax, log10(l23ShVols), 'BinEdges', curBinEdges);

hists = flip(ax.Children);
set(hists, 'Normalization', 'probability');

xlabel(ax, 'Spine head volume [log10(µm³)]');
ylabel(ax, 'Probability');

leg = legend(hists, { ...
    sprintf('L4 (n = %d)', numel(l4ShVols)), ...
    sprintf('L2/3 (n = %d)', numel(l23ShVols))});
leg.Location = 'EastOutside';

connectEM.Figure.config(fig, info);
fig.Position(3:4) = [460, 225];

%% Plot ASI areas
clear cur*;
curBinEdges = linspace(-2, 0.5, 26);

fig = figure();
ax = axes(fig);

hold(ax, 'on');
histogram(ax, log10(l4AsiT.area), 'BinEdges', curBinEdges);
histogram(ax, log10(l23AsiT.area), 'BinEdges', curBinEdges);

hists = flip(ax.Children);
set(hists, 'Normalization', 'probability');

xlabel(ax, 'ASI area [log10(µm²)]');
ylabel(ax, 'Probability');

leg = legend(hists, { ...
    sprintf('L4 (n = %d)', height(l4AsiT)), ...
    sprintf('L2/3 (n = %d)', height(l23AsiT))});
leg.Location = 'EastOutside';

connectEM.Figure.config(fig, info);
fig.Position(3:4) = [460, 225];

%% Correlate spine head volume with ASI area
clear cur*;
curLimX = [-4, +1];
curLimY = [-3, 1];

fig = figure();
subplot(1, 2, 1);
scatter(log10(l4ShVols(l4AsiT.shId)), log10(l4AsiT.area), '.');
curFit = fitlm(log10(l4ShVols(l4AsiT.shId)), log10(l4AsiT.area)) %#ok
hold('on'); plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');
axis('square'); xlim(curLimX); ylim(curLimY);
xlabel('Spine head volume [log10(µm³)]');
ylabel('ASI area [log10(µm²)]');

subplot(1, 2, 2);
scatter(log10(l23ShVols(l23AsiT.shId)), log10(l23AsiT.area), '.');
curFit = fitlm(log10(l23ShVols(l23AsiT.shId)), log10(l23AsiT.area)) %#ok
hold('on'); plot(curLimX(:), curFit.predict(curLimX(:)), 'k--');
curLeg = legend(sprintf('y = %.2f + %.2fx (R² = %.2f)', ...
    curFit.Coefficients.Estimate, curFit.Rsquared.Ordinary));
set(curLeg, 'Box', 'Off', 'Location', 'South');
axis('square'); xlim(curLimX); ylim(curLimY);

connectEM.Figure.config(fig, info);
fig.Position(3:4) = [555, 250];

%% Plot spine head volumes (attached vs. unattached)
clear cur*;
curBinEdges = linspace(-3, 0, 41);

fig = figure();
ax = axes(fig);

hold(ax, 'on');
histogram(ax, log10(l23ShVols( ...
    l23ShAttached > 0)), 'BinEdges', curBinEdges);
histogram(ax, log10(l23ShVols( ...
    l23ShAttached == 0)), 'BinEdges', curBinEdges);

hists = flip(ax.Children);
set(hists, 'Normalization', 'probability');

xlabel(ax, 'Spine head volume [log10(µm³)]');
ylabel(ax, 'Probability');

leg = legend(hists, { ...
    sprintf('L2/3 attached (n = %d)', sum(l23ShAttached > 0)), ...
    sprintf('L2/3 unattached (n = %d)', sum(l23ShAttached == 0))});
leg.Location = 'EastOutside';

connectEM.Figure.config(fig, info);
fig.Position(3:4) = [545, 225];
