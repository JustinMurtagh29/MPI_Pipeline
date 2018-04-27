% Classifies dendrites into
% * smooth dendrites (based on spine density)
% * apical dendrites (based on manual annotations)
% * other dendrites
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
adFile = '/tmpscratch/sahilloo/L4/dataPostSyn/dendritesADState_v2.mat';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

[outDir, outFile] = fileparts(dendFile);
outFile = fullfile(outDir, sprintf('%s_classified.mat', outFile));
clear outDir;

% Set path to export NML file with conflicts
confNmlFile = '';

% NML file resolving conflicts
annNmlFile = fullfile( ...
    fileparts(mfilename('fullpath')), 'annotations', ...
    '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__b3a207.nml');

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex
maxSpinesPerUm = 0.4;

info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Apical dendrites
adAgglos = load(adFile);
adAgglos = adAgglos.agglos;

% Dendrite trunks (i.e., prior to spine attachment)
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Spine heads
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Dendrites
dend = load(dendFile);

dendIds = dend.indAIS | dend.indSomata | dend.indWholeCells;
dendIds = find(dend.indBigDends & ~dendIds);

dendAgglos = dend.dendrites(dendIds);
dendAgglos = Agglo.fromSuperAgglo(dendAgglos);

%% Smooth dendrites
Util.log('Calculating spine density');

trunkLens = ...
    connectEM.Dendrite.calculatePathLengths(param, dendAgglos, trunks);
spineCounts = ...
    connectEM.Dendrite.calculateSpineCount(param, dendAgglos, shAgglos);
spineDensity = spineCounts ./ (trunkLens / 1E3);

sdIds = dendIds(spineDensity < maxSpinesPerUm);

%% Apical dendrites
maxSegId = Seg.Global.getMaxSegId(param);
dendLUT = Agglo.buildLUT(maxSegId, dendAgglos, dendIds);

adIds = cellfun(@(ids) mode(nonzeros(dendLUT(ids))), adAgglos);
adIds = unique(adIds(~isnan(adIds)));

%% Export conflicts
confIds = intersect(sdIds, adIds);
sdIds = setdiff(sdIds, confIds);
adIds = setdiff(adIds, confIds);

if ~isempty(confNmlFile)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    for curId = reshape(confIds, 1, [])
        curAgglo = dend.dendrites(curId);
        skel = skel.addTree( ...
            sprintf('Dendrite %d', curId), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
    end

    skel.write(confNmlFile);
end

%% Solve conflicts using MH's annotations
nml = slurpNml(annNmlFile);
trees = NML.buildTreeTable(nml);

% Find agglomerate IDs and get rid of invalid trees
trees.dendId = regexpi(trees.name, '^Dendrite\W+(\d+)', 'tokens', 'once');
trees(cellfun(@isempty, trees.dendId), :) = [];
trees.dendId = cellfun(@str2double, trees.dendId);

% Sanity check
assert(all(ismember(trees.dendId, dendIds)));

trees.targetClass(:) = {'OtherDendrite'};
trees.targetClass(contains( ...
    trees.name, 'mhSD', 'IgnoreCase', true)) = {'SmoothDendrite'};
trees.targetClass(contains( ...
    trees.name, 'mhAD', 'IgnoreCase', true)) = {'ApicalDendrite'};
trees.targetClass = categorical(trees.targetClass);

sdIds = union(sdIds, trees.dendId(trees.targetClass == 'SmoothDendrite'));
adIds = union(adIds, trees.dendId(trees.targetClass == 'ApicalDendrite'));

% Sanity check
assert(isempty(intersect(sdIds, adIds)));

%% Build output
out = dend;
out.indSmoothies = false(size(out.dendrites));
out.indSmoothies(sdIds) = true;
out.indApicals = false(size(out.dendrites));
out.indApicals(adIds) = true;

% Build `targetClass` categorical
out.targetClass = cell(size(out.dendrites));
out.targetClass(:) = {'Ignore'};
out.targetClass(out.indBigDends) = {'OtherDendrite'};
out.targetClass(out.indSomata) = {'Somata'};
out.targetClass(out.indWholeCells) = {'WholeCell'};
out.targetClass(out.indAIS) = {'AxonInitialSegment'};
out.targetClass(out.indApicals) = {'ApicalDendrite'};
out.targetClass(out.indSmoothies) = {'SmoothDendrite'};
out.targetClass = categorical(out.targetClass);

out = rmfield(out, 'info');
out = orderfields(out);
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
