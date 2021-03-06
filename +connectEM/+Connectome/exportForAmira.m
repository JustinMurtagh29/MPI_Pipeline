% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
outputDir = '/tmpscratch/amotta/l4/2019-09-27-connectome-isosurfaces';

isoParams = { ...
    'reduce', 0.05, ...
    'smoothSizeHalf', 4, ...
    'smoothWidth', 8};

segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.seg = segParam;

maxSegId = Seg.Global.getMaxSegId(param);
conn = connectEM.Connectome.load(param, connFile);

dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites;

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Connectome that was used for consistency analysis
% See commit 34dee843596df4ad59a68bb2b8cf0d03edf11038
% of +connectEM/+Consistency/+Script/buildAxonSpineInterfaces.m
linConn = connectEM.Consistency.loadConnectomePaper(param);
linConn = connectEM.Connectome.prepareForSpecificityAnalysis(linConn);

%% Generate isosurfaces of axons
% See also +connectEM/+Axon/+Script/calculateIsosurfaces.m
clear cur*;

% Remove soma segments from axons
curSomaMask = conn.dendrites(conn.denMeta.targetClass == 'Somata');
curSomaMask = logical(Agglo.buildLUT(maxSegId, curSomaMask));

curAxons = cellfun( ...
    @(segIds) segIds(~curSomaMask(segIds)), ...
    conn.axons, 'UniformOutput', false);

curOutDir = fullfile(outputDir, 'axons');
assert(not(exist(curOutDir, 'dir')));
assert(mkdir(curOutDir));

Visualization.exportAggloToAmira( ...
    param, curAxons, curOutDir, isoParams{:});

curInfoFile = fullfile(curOutDir, 'info.mat');
Util.save(curInfoFile, info);
Util.protect(curOutDir);

%% Generate isosurfaces of dendrites
clear cur*;

curOutDir = fullfile(outputDir, 'dendrites');
assert(not(exist(curOutDir, 'dir')));
assert(mkdir(curOutDir));

curDendrites = conn.dendrites;
Visualization.exportAggloToAmira( ...
    param, curDendrites, curOutDir, isoParams{:});

curInfoFile = fullfile(curOutDir, 'info.mat');
Util.save(curInfoFile, info);
Util.protect(curOutDir);

%% Generate isosurfaces of linearized axons used for consistency analysis
% See +connectEM/+Connectome/plotSynapseSizeConsistency.m
% commit 25b3a8f663b86656cf9b149b63fd50e84d68fa81, and
% +connectEM/+Consistency/+Script/buildAxonSpineInterfaces.m
% commit 34dee843596df4ad59a68bb2b8cf0d03edf11038
clear cur*;

% Remove soma segments from axons
curSomaMask = linConn.dendrites(linConn.denMeta.targetClass == 'Somata');
curSomaMask = logical(Agglo.buildLUT(maxSegId, curSomaMask));

curAxons = cellfun( ...
    @(segIds) segIds(~curSomaMask(segIds)), ...
    linConn.axons, 'UniformOutput', false);

curOutDir = fullfile(outputDir, 'axons-linear');
assert(not(exist(curOutDir, 'dir')));
assert(mkdir(curOutDir));

Visualization.exportAggloToAmira( ...
    param, curAxons, curOutDir, isoParams{:});

curInfoFile = fullfile(curOutDir, 'info.mat');
Util.save(curInfoFile, info);
Util.protect(curOutDir);

%% Generate NML file with skeletons of manual spine head attachments
clear cur*;

curShLUT = Agglo.buildLUT(maxSegId, shAgglos);

curSkels = dendrites(conn.denMeta.parentId);
curDendIds = reshape(1:numel(curSkels), [], 1);
[~, curSkels] = Superagglos.splitIntoAgglosAndFlights(curSkels);

curMask = not(cellfun(@isempty, curSkels));
curDendIds = curDendIds(curMask);
curSkels = curSkels(curMask);

curCounts = repelem( ...
    1:numel(curSkels), ...
    cellfun(@numel, curSkels));
curDendIds = curDendIds(curCounts);
curSkels = reshape(cat(2, curSkels{:}), [], 1);

curShIds = arrayfun( ...
    @(s) s.nodes(s.nodes(:, 4) > 0, 4), ...
    curSkels, 'UniformOutput', false);
curShIds = cellfun( ...
    @(s) reshape(setdiff(curShLUT(s), 0), [], 1), ...
    curShIds, 'UniformOutput', false);

curCounts = repelem( ...
    1:numel(curShIds), ...
    cellfun(@numel, curShIds));

curDendIds = curDendIds(curCounts);
curSkels = curSkels(curCounts);
curShIds = cat(1, curShIds{:});

out = struct;
out.info = info;
out.skels = curSkels(:);
out.shIds = curShIds(:);
out.dendIds = curDendIds(:);
