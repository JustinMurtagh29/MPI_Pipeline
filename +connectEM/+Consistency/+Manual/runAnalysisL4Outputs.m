% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

curNmlDir = fullfile(fileparts(mfilename('fullpath')), 'annotations');
adNmlFile = fullfile(curNmlDir, 'l4-ad-spine-synapses.nml');
pdNmlFile = fullfile(curNmlDir, 'l4-pd-spine-synapses.nml');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
clear cur*;

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
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Show results
clear cur*;

outT = table;
[outT.targetClass, ~, curIds] = unique(asiT.targetClass);
outT.areas = accumarray(curIds, asiT.area, [], @(areas) {areas(:)});
outT.synCount = cellfun(@numel, outT.areas);
outT.meanArea = cellfun(@mean, outT.areas);
outT.medianArea = cellfun(@median, outT.areas);

fprintf('Summary\n\n');
disp(outT);

% Statistical tests
statT = struct;
statT(1).name = 'Mann-Whitney-Wilcoxon rank sum test';
statT(1).pVal = ranksum(outT.areas{:});

statT(2).name = 'Kolmogorov-Smirnov test';
[~, statT(2).pVal] = kstest2(outT.areas{:});

curData = cellfun(@log10, outT.areas, 'UniformOutput', false);
statT(3).name = 'Student''s t-test in log10';
[~, statT(3).pVal] = ttest2(curData{:});

fprintf('Statistical tests\n\n');
statT = struct2table(statT, 'AsArray', true);
disp(statT);

%% Histograms
clear cur*;

curConfigs = struct;
curConfigs(1).xlabel = 'log_{10}(ASI area [µm²])';
curConfigs(1).binEdges = linspace(-1.5, 0.5, 21);
curConfigs(1).data = cellfun(@log10, outT.areas, 'UniformOutput', false);

curConfigs(2).xlabel = 'ASI area [µm²]';
curConfigs(2).binEdges = linspace(0, 0.8, 17);
curConfigs(2).data = outT.areas;

for curConfig = curConfigs
    curFig = figure();
    curFig.Position(3:4) = [290, 290];
    
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');

    cellfun( ...
        @(values) histogram( ...
            curAx, values, ...
            'BinEdges', curConfig.binEdges), ...
        curConfig.data);

    xlabel(curAx, curConfig.xlabel);
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
end
