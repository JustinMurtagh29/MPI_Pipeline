% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% TODO(amotta): FIX FIX FIX!
nmlDir = '/mnt/mpibr/data/Personal/mottaa/L4/2019-02-19-Proofreading-L4-Output-Synapses';
adNmlFile = fullfile(nmlDir, 'l4-onto-ad-spine-synapses_merger-mode_v1.nml');
pdNmlFile = fullfile(nmlDir, 'l4-onto-pd-spine-synapses_merger-mode_v1.nml');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToCentroidMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

%% Load NML files
curLoad = @connectEM.Consistency.Manual.loadAnnotations;
curUnpack = @(v) v{1};

curAdAnnT = curUnpack(curLoad(param, adNmlFile));
curAdAnnT.targetClass(:) = categorical({'ApicalDendrite'});

curPdAnnT = curUnpack(curLoad(param, pdNmlFile));
curPdAnnT.targetClass(:) = categorical({'ProximalDendrite'});

annT = vertcat(curAdAnnT, curPdAnnT);

%% Calculate ASI areas
clear cur*;
curAxons = annT.preSegIds;
curSpineHeads = annT.postSegIds;

Agglo.check(vertcat(curAxons, curSpineHeads));

asiT = annT(:, {'id', 'targetClass'});
asiT.preAggloId(:) = 1:numel(curAxons);
asiT.postAggloId(:) = nan;
asiT.shId(:) = 1:numel(curSpineHeads);

curWmean = @(v, w) sum(v .* (w / sum(w)), 1);
asiT.pos = cell2mat(cellfun( ...
    @(ids) curWmean(segPoints(ids, :), segSizes(ids, :)), ...
    curSpineHeads, 'UniformOutput', false));
asiT.pos = round(asiT.pos);

asiT.area = ...
    connectEM.Consistency.buildAxonSpineInterfaceAreas( ...
        param, curAxons, curSpineHeads, asiT);

%% Show results
clear cur*;
curBinEdges = linspace(-1.5, 0.5, 21);

outT = table;
[outT.targetClass, ~, curIds] = unique(asiT.targetClass);
outT.synCount = accumarray(curIds, 1);
outT.meanArea = accumarray(curIds, asiT.area, [], @mean);
outT.medianArea = accumarray(curIds, asiT.area, [], @median);

disp(outT);

% Histograms
curFig = figure();
curFig.Position(3:4) = [290, 290];
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:height(outT)
    curAreas = asiT.area(curIds == curIdx);
    histogram(curAx, log10(curAreas), 'BinEdges', curBinEdges);
end

xlabel(curAx, 'log_{10}(ASI area [µm²])');
ylabel(curAx, 'Occurences');

curLeg = legend(curAx, ...
    arrayfun( ...
        @(t, n) sprintf('%s (n = %d)', t, n), ...
        outT.targetClass, outT.synCount, ...
        'UniformOutput', false), ...
    'Location', 'SouthOutside');

connectEM.Figure.config(curFig, info);

curHists = findobj(curFig, 'type', 'histogram');
set(curHists, 'DisplayStyle', 'bar', 'EdgeColor', 'none');
