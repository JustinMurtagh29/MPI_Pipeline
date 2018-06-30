% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

beforeFile = fullfile(rootDir, 'aggloState', 'dendrites_03_v2.mat');
afterFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

before = load(beforeFile);
after = load(afterFile);

%% Prepare look-up tables
% TODO(amotta): With or without `indBig` filter?
beforeAgglos = Agglo.fromSuperAgglo(before.dendrites);
afterAgglos = Agglo.fromSuperAgglo(after.wholeCells);

beforeLUT = Agglo.buildLUT(maxSegId, beforeAgglos);
afterLUT = Agglo.buildLUT(maxSegId, afterAgglos);

%% Estimate number of solved splits
beforeAgglosCollected = cellfun( ...
    @(segIds) setdiff(beforeLUT(segIds), 0), ...
    afterAgglos, 'UniformOutput', false);
numSplitsPerCell = cellfun(@numel, beforeAgglosCollected);
numSplits = sum(numSplitsPerCell) %#ok

%% Estimate number of solved mergers
beforeAgglosModified = cellfun(@(ids) arrayfun(@(id) ...
        ~isscalar(unique(afterLUT(beforeAgglos{id}))), ids), ...
	beforeAgglosCollected, 'UniformOutput', false);
numMergersPerCell = cellfun(@sum, beforeAgglosModified);
numMergers = sum(numMergersPerCell) %#ok
