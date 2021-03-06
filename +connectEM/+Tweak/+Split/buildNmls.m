% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aggloFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');
outDir = '/home/amotta/Desktop/axon-nmls';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

agglos = load(aggloFile);
agglos = agglos.axons(agglos.indBigAxons);

%% Select agglomerates to export
segCount = arrayfun(@(a) sum(not(isnan(a.nodes(:, 4)))), agglos);
[segCount, aggloIds] = sort(segCount, 'descend');

% Ignore the blood vessel
segCount(1) = [];
aggloIds(1) = [];

aggloIds = aggloIds(51:100);

%% Generate NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDescription = sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash);
numDigits = ceil(log10(1 + numel(aggloIds)));

for curIdx = 1:numel(aggloIds)
    curId = aggloIds(curIdx);
    curAgglo = agglos(curId);
    
    curSkel = connectEM.Tweak.buildNml(curAgglo, curId, skel);
    curSkelName = sprintf('%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outDir, curSkelName));
end
