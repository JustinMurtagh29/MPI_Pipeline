% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'spine_heads_and_attachment_03.mat');

% position of spine heads
shCandPos = 1 + [ ...
    765, 2231, 457;
    1944, 8015, 2821;
    964, 1256, 2405;
    2973, 5692, 2566];

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

maxSegId = Seg.Global.getMaxSegId(param);
shCandSegIds = Seg.Global.getSegIds(param, shCandPos);

%% try to find candidates in spine head agglomerates
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
shCandAggloIds = shLUT(shCandSegIds)
