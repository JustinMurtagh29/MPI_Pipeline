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
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03.mat');
adFile = '/tmpscratch/sahilloo/L4/dataPostSyn/dendritesADState_v2.mat';

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

dend = load(dendFile);

adAgglos = load(adFile);
adAgglos = adAgglos.agglos;

%% Build dendrite look-up table
maxSegId = Seg.Global.getMaxSegId(param);

dendIds = dend.indAIS | dend.indWholeCells | dend.indSomata;
dendIds = find(dend.indBigDends & ~dendIds);

dendAgglos = Agglo.fromSuperAgglo(dend.dendrites(dendIds));
dendLUT = Agglo.buildLUT(maxSegId, dendAgglos, dendIds);

%% Find apical dendrites
adDendIds = cellfun(@(ids) mode(nonzeros(dendLUT(ids))), adAgglos);
adDendIds = unique(adDendIds(~isnan(adDendIds)));

%% TODO(amotta): Find smooth dendrites

%% TODO(amotta): Save results
