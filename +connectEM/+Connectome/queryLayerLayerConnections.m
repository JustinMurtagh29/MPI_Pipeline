% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';
evalNmlDir = connectEM.Consistency.Manual.getDir(fullfile('annotations', 'l4-spine-synapses'));

l4ConnRunId = '20190221T112510';
outDir = '/home/amotta/Desktop';

checkedSynIds = [ ...
    ... All primary spine synapses from L4 → AD (19.02.2019)
    7975, 68992, 77903, 185118, 186018, 97233, 119754, ...
    122684, 271912, 300946, 344896, 294170, 140701, 202773, ...
    ... Randomly selected primary spine synapses from L4 → L4 (20.02.2019)
    97662, 109112, 72299, 101803, 228602, 17999, ...
    275717, 163160, 44870, 163175, 33945, 94428];

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

curData = load(outputMapFile);
axonData = curData.axonData;

conn = curData.info.param.connFile;
conn = connectEM.Connectome.load(param, conn);

[curDir, curFile] = fileparts(outputMapFile);
curAsiFile = sprintf('%s__%s_connectome.mat', curFile, l4ConnRunId);
curData = load(fullfile(curDir, curAsiFile));

l4SynT = curData.synT;
l4AsiT = curData.asiT;
l4AsiT = l4AsiT(l4AsiT.area > 0, :);

%% Prepare synapse table
l4SynT.targetClass = conn.denMeta.targetClass(l4SynT.postAggloId);
[~, curIds] = ismember(l4AsiT.id, l4SynT.id);

l4SynT.type(:) = categorical({'Shaft'});
l4SynT.type(curIds) = l4AsiT.type;

%% Select random synapses to query
clear cur*;
rng(0);

curSynT = l4SynT;
curSynT = curSynT(curSynT.type == 'PrimarySpine', :);

% 20.02.2019: All L4 → AD and 20 L4 → WC synapses
querySynIds = checkedSynIds;
curSynT(ismember(curSynT.id, querySynIds), :) = [];

% 22.02.2019: Select 40 more L4 → WC synapses
curSynIds = curSynT.id(curSynT.targetClass == 'WholeCell');
curSynIds = curSynIds(randperm(numel(curSynIds)));
curSynIds = curSynIds(1:40);

querySynIds = union(querySynIds, curSynIds, 'stable');
curSynT(ismember(curSynT.id, querySynIds), :) = [];

% 27.02.2019: Select 40 L4 → OD synapses
curSynIds = curSynT.id(curSynT.targetClass == 'OtherDendrite');
curSynIds = curSynIds(randperm(numel(curSynIds)));
curSynIds = curSynIds(1:40);

querySynIds = curSynIds;
curSynT(ismember(curSynT.id, querySynIds), :) = [];

%% Generate table for Heiko
clear cur*;

queryT = table;
queryT.synId = querySynIds(:);

[~, curIds] = ismember(queryT.synId, l4SynT.id);
[~, queryT.nmlFile] = arrayfun( ...
    @(axonId) fileparts(axonData(axonId).nmlFile), ...
    l4SynT.preAggloId(curIds), 'UniformOutput', false);

format long;
fprintf('Queries\n\n');
disp(queryT);

%% Evaluate queries
clear cur*;

curCalcAsiArea = @connectEM.Consistency.Calibration.calculateAsiArea;

evalNmlFiles = dir(fullfile(evalNmlDir, '*.nml'));
evalNmlFiles(cat(1, evalNmlFiles.isdir)) = [];
evalNmlFiles = fullfile(evalNmlDir, {evalNmlFiles.name});

evalAsiT = [];

for curIdx = 1:numel(evalNmlFiles)
    curEvalNmlFile = evalNmlFiles{curIdx};
   [~, curEvalNmlName] = fileparts(curEvalNmlFile);
   
    curNml = slurpNml(curEvalNmlFile);
    curTrees = NML.buildTreeTable(curNml);
    curGroups = struct2table(curNml.groups, 'AsArray', false);
    
    % Find groups with synapse annotations
    curAsiT = table;
    curAsiT.groupId = curGroups.id;
    curAsiT.name = curGroups.name;
    
    curAsiT.id = regexpi( ...
        curGroups.name, '^synapse (\d+)', 'tokens', 'once');
    curAsiT = curAsiT(cellfun(@isscalar, curAsiT.id), :);
    curAsiT.id = cellfun(@(id) str2double(id{1}), curAsiT.id);
    
    curAsiT.isCorrect = not(contains( ...
        curAsiT.name, 'wrong', 'IgnoreCase', true));
    
    % Measure synapse areas
    curAsiT.area(:) = nan;
    
    for curSynIdx = 1:height(curAsiT)
        if not(curAsiT.isCorrect(curSynIdx)); continue; end
        
        curGroupId = curAsiT.groupId(curSynIdx);
        curSynTrees = curTrees(curTrees.groupId == curGroupId, :);
        
        try
            curArea = curCalcAsiArea(param.raw.voxelSize, curSynTrees);
        catch
            % NOTE(amotta): The ASI area calculation routine throws an
            % error in case the tracings are invalid (e.g., because a slice
            % has been skipped). Let's ignore these synapses and leave
            % their area at `nan` for now...
            warning( ...
                'Group "%s" of in %s.nml is invalid', ...
                curAsiT.name{curSynIdx}, curEvalNmlName);
            continue;
        end
        
        curAsiT.area(curSynIdx) = curArea;
    end
    
    % Clean-up
    curAsiT.name = [];
    curAsiT.groupId = [];
    
    % HACK(amotta): Lazy initialization
    if isempty(evalAsiT); evalAsiT = curAsiT([], :); end
    evalAsiT = cat(1, evalAsiT, curAsiT);
end

[~, evalAsiT.targetClass] = ismember(evalAsiT.id, l4SynT.id);
evalAsiT.targetClass = l4SynT.targetClass(evalAsiT.targetClass);

%% Quantitative evaluation
clear cur*;

curOutT = table;
[curOutT.targetClass, ~, curIds] = ...
    unique(evalAsiT.targetClass);
curOutT.areas = accumarray( ...
    curIds, evalAsiT.area, [], ...
    @(areas) {rmmissing(areas(:))});

curOutT.synCount = cellfun(@numel, curOutT.areas);
curOutT.meanArea = cellfun(@mean, curOutT.areas);
curOutT.medianArea = cellfun(@median, curOutT.areas);

fprintf('Summary\n\n');
disp(curOutT);

% Statistical tests
curStatT = struct;
curStatT(1).name = 'Mann-Whitney-Wilcoxon rank sum test';
curStatT(1).pVal = ranksum(curOutT.areas{:});

curStatT(2).name = 'Kolmogorov-Smirnov test';
[~, curStatT(2).pVal] = kstest2(curOutT.areas{:});

curData = cellfun(@log10, curOutT.areas, 'UniformOutput', false);
curStatT(3).name = 'Student''s t-test in log10';
[~, curStatT(3).pVal] = ttest2(curData{:});

fprintf('Statistical tests\n\n');
curStatT = struct2table(curStatT, 'AsArray', true);
disp(curStatT);

%% Plots for evaluation
clear cur*;

curConfigs = struct;
curConfigs(1).xlabel = 'log_{10}(ASI area [µm²])';
curConfigs(1).binEdges = linspace(-1.5, 0.5, 21);
curConfigs(1).transform = @log10;

curConfigs(2).xlabel = 'ASI area [µm²]';
curConfigs(2).binEdges = linspace(0, 0.8, 17);
curConfigs(2).transform = @(id) id;

curClasses = unique(evalAsiT.targetClass);

for curConfig = curConfigs
    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');

    for curClass = reshape(curClasses, 1, [])
        curAreas = evalAsiT.area(evalAsiT.targetClass == curClass);
        curAreas = curConfig.transform(curAreas);
        
        histogram( ...
            curAx, curAreas, ...
            'BinEdges', curConfig.binEdges, ...
            'Normalization', 'probability');
    end
    
    xlabel(curAx, curConfig.xlabel);
    ylabel(curAx, 'Probability');
    
    curHists = flip(findobj(curFig, 'type', 'histogram'));
    
    curLeg = arrayfun( ...
        @(c, h) sprintf('%s (n = %d)', c, numel(h.Data)), ...
        curClasses, curHists, 'UniformOutput', false);
    curLeg = legend(curAx, curLeg, 'Location', 'SouthOutside');
    
    curFig.Position(3:4) = [220, 260];
    connectEM.Figure.config(curFig, info);
    
    set(curHists, 'DisplayStyle', 'bar', 'EdgeColor', 'none');
end
