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
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_autoSpines_v1.mat');

[outDir, outFile] = fileparts(dendFile);
outFile = fullfile(outDir, sprintf('%s_classified_v2.mat', outFile));
clear outDir;

% Set path to export NML file with conflicts
confNmlFile = '';

% NML file resolving conflicts
annNmlFile = fullfile( ...
    fileparts(mfilename('fullpath')), 'annotations', ...
    'sd-ad-conflict-resolution_auto-pre-robo.nml');

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex

% NOTE(amotta): This automated pre-RoboEM state only contains spines that
% got attached by the greedy walk-based spine attachment method. This
% systematically lowers the spine density compared to the dendrite
% reconstruction that also contained manually attached spines.
%   To prevent a systematic bias towards smooth dendrites, let's lower the
% spine density threshold that is used to define smooth dendrites. We do
% this by multiplying the threshold with the ration of automatically to
% automatically-or-manually attached spines.

% NOTE(amotta): These numbers were obtained by
% connectEM.Number.spineAttachment
% git@gitlab.mpcdf.mpg.de:connectomics/pipeline.git f33c7c2f786401135c7f3715f802beb519c2b56b
% amotta@m-01522. MATLAB 9.3.0.713579 (R2017b). 02-Jun-2020 14:17:05
numAutoAttachedSpines = 229990; % sum(shT.autoAttached)
numAutoOrManuallyAttachedSpines = 327958; % sum(shT.attached)
maxSpinesPerUm = 0.4 * (numAutoAttachedSpines / numAutoOrManuallyAttachedSpines);

info = Util.runInfo();
Util.showRunInfo(info);

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
    skel = Skeleton.setDescriptionFromRunInfo(skel, info);

    for curId = reshape(confIds, 1, [])
        curAgglo = dend.dendrites(curId);
        skel = skel.addTree( ...
            sprintf('Dendrite %d', curId), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
    end

    skel.write(confNmlFile);
end

%% Solve conflicts using MH's annotations
if ~isempty(annNmlFile)
    nml = slurpNml(annNmlFile);
    trees = NML.buildTreeTable(nml);

    % Find agglomerate IDs and get rid of invalid trees
    trees.dendId = regexpi(trees.name, ...
        '^Dendrite\W+(\d+)', 'tokens', 'once');
    trees(cellfun(@isempty, trees.dendId), :) = [];
    trees.dendId = cellfun(@str2double, trees.dendId);
    
    curSdMask = contains(trees.name, {'mhSD', 'amSD'}, 'IgnoreCase', true);
    curAdMask = contains(trees.name, {'mhAD', 'amAD'}, 'IgnoreCase', true);

    % Sanity check
    assert(all(ismember(trees.dendId, dendIds)));

    trees.targetClass(:) = {'OtherDendrite'};
    trees.targetClass(curSdMask) = {'SmoothDendrite'};
    trees.targetClass(curAdMask) = {'ApicalDendrite'};
    trees.targetClass = categorical(trees.targetClass);

    sdIds = union(sdIds, trees.dendId( ...
        trees.targetClass == 'SmoothDendrite'));
    adIds = union(adIds, trees.dendId( ...
        trees.targetClass == 'ApicalDendrite'));
else
    error('Require NML file for conflict resolution to continue');
end

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
