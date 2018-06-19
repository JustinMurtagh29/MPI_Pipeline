% This script generates NML files of soma-alignmed whole cells. This code
% was originally part of `connectEM.Connectome.plotWholeCellInputs`.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/soma-aligned-cells';

wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segMass = Seg.Global.getSegToSizeMap(param);
segPoint = Seg.Global.getSegToPointMap(param);

wcData = load(wcFile);
somaData = load(somaFile);

%% Complete whole cells
wcT = table;
wcT.id = wcData.idxWholeCells(wcData.indWholeCells);
wcT.agglo = wcData.dendrites(wcData.indWholeCells);
SuperAgglo.check(wcT.agglo);

% Find corresponding somata
[~, somaIds] = ismember(wcT.id, somaData.idxSomata);
wcT.somaAgglo = somaData.dendrites(somaIds);

calcSomaPos = @(n) ...
    sum(segMass(n(:, 4)) .* n(:, 1:3), 1) ....
    ./ sum(segMass(n(:, 4)));
wcT.somaPos = cell2mat(arrayfun( ...
    @(a) calcSomaPos(a.nodes), ...
    wcT.somaAgglo, 'UniformOutput', false));

% NOTE(amotta): There is a bug in the soma super-agglomerates which allows
% them to be disconnected. Let's fix this by introducing random edges. This
% is not a problem since we do not care about distances within the soma.
wcT.somaAgglo = SuperAgglo.connect(wcT.somaAgglo);
SuperAgglo.check(wcT.somaAgglo);

% Merge somata with whole cells
wcT.agglo = arrayfun(@SuperAgglo.merge, wcT.agglo, wcT.somaAgglo);
wcT.agglo = SuperAgglo.clean(wcT.agglo);

% Sort by increasing node count
wcT.nodeCount = arrayfun(@(a) size(a.nodes, 1), wcT.agglo);
wcT = sortrows(wcT, 'nodeCount', 'ascend');
wcT.nodeCount = [];

%% Build NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

wcCount = height(wcT);
numDigits = ceil(log10(1 + wcCount));

for curIdx = 1:wcCount
    curWc = table2struct(wcT(curIdx, :));
    
    curAgglo = curWc.agglo;
    curAgglo.nodes(:, 1:3) = ...
        curAgglo.nodes(:, 1:3) ...
      - curWc.somaPos;
    
    curSkel = Superagglos.toSkel(curAgglo, skel);
    curSkel.names{1} = sprintf(...
        '%0*d. Whole cell %d', numDigits, curIdx, curWc.id);
    
    curSkelName = sprintf( ...
        '%0*d_whole-cell-%d.nml', numDigits, curIdx, curWc.id);
    curSkel.write(fullfile(outputDir, curSkelName));
end
