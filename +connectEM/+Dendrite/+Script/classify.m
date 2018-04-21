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
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2_auto.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03.mat');

nmlFile = '/home/amotta/Desktop/sd-ad-conflicts.nml';

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

%% TODO(amotta): Build output
% Sanity check
% assert(isempty(intersect(sdIds, adIds)));

%% Export conflicts
if ~isempty(nmlFile)
    confIds = intersect(sdIds, adIds);

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

    skel.write(nmlFile);
end
