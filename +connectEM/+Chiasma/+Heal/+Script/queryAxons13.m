% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_13_a.mat');

outputDir = fullfile( ...
    rootDir, 'chiasmaSplitHealing', ...
    '20171207T112927-healing-axons-13');
assert(logical(exist(outputDir, 'dir')));

% segmentation for dynamic stopping of queries
segRoot = '/tmpscratch/amotta/l4/2017-12-07-segmentation-for-axons-13';

info = Util.runInfo();

%% loading parameters
paramFile = fullfile(rootDir, 'allParameter.mat');
param = load(paramFile, 'p');
param = param.p;

%% loading axons
axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% generate task and webKNOSSOS task definitions
endings = connectEM.Chiasma.Heal.buildEndings(axons);
taskDefFile = fullfile(outputDir, sprintf('%s_flightTasks.txt', runId));

[taskDefs, endings] = ...
    connectEM.Chiasma.Heal.generateTasks( ...
        param, endings, axons, taskDefFile);

%%
out = struct;
out.info = info;
out.endings = endings;
out.taskDefs = taskDefs;
out.taskDefFile = taskDefFile;

outFile = sprintf('%s_taskGeneration.mat', runId);
Util.saveStruct(fullfile(outputDir, outFile), out);

%% generate segmentation for dynamic stopping
%{
maxSegId = Seg.Global.getMaxSegId(param);
agglos = Superagglos.getSegIds(axons);
mapping = Agglo.buildLUT(maxSegId, agglos);

seg = struct;
seg.root = fullfile(segRoot, '1');
seg.prefix = param.seg.prefix;
seg.backend = 'wkwrap';

% Initialize WKW dataset
wkwInit('new', seg.root, 32, 32, 'uint32', 1);
Seg.Global.applyMappingToSegmentation(param, mapping, seg);
compressSegmentation(seg);
%}