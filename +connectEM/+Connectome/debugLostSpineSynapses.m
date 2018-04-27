% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');

connFiles = { ...
    'connectome_axons_18_a_ax_spine_syn_clust.mat';
    'connectome_axons-19-a_dendrites-wholeCells-03-classified_spine-syn-clust.mat'};
connFiles = fullfile(rootDir, 'connectomeState', connFiles);

shFiles = { ...
    'dendrites_wholeCells_01_auto.mat';
    'dendrites_wholeCells_02_v2_auto.mat'};
shFiles = fullfile(rootDir, 'aggloState', shFiles);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

syn = load(synFile);
connOld = load(connFiles{1});
connNew = load(connFiles{2});

shOld = load(shFiles{1});
shNew = load(shFiles{2});

%% Spine heads + attachment
% Compare spine head agglomerates
assert(isequal(shOld.shAgglos, shNew.shAgglos));

oldLUT = Agglo.buildLUT(maxSegId, shOld.dendAgglos);
newLUT = Agglo.buildLUT(maxSegId, shNew.dendAgglos);

% Compare attachment
oldMask = (shOld.attached ~= 0);
oldMask2 = cellfun(@(ids) any(oldLUT(ids)), shOld.shAgglos);

sum(oldMask & ~oldMask2)
sum(oldMask2 & ~oldMask)

newMask = (shNew.attached ~= 0);
newMask2 = cellfun(@(ids) any(newLUT(ids)), shNew.shAgglos);

sum(newMask & ~newMask2)
sum(newMask2 & ~newMask)

sum(oldMask & ~newMask)
sum(newMask & ~oldMask)
