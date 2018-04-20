% Make sure that axons_18_a and axons_18_b are indeed equivalent in terms
% of the segments that they contain.
%
% Please note that axons 18a and 18b are not *identical*. The order in
% which super-agglomerates are merged (and redundant parts are removed) was
% improved from 18a to 18b. This change in order means that that different
% nodes were identified as redundant and removed.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aFile = fullfile(rootDir, 'aggloState', 'axons_18_a.mat');
bFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');

% Set to 
nmlId = 281;
nmlDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

a = load(aFile);
b = load(bFile);

%% Checking segment equivalence classes on their own
aAgglos = Agglo.fromSuperAgglo(a.axons);
Agglo.check(aAgglos);

bAgglos = Agglo.fromSuperAgglo(b.axons);
Agglo.check(aAgglos);

%% Check equivalence of segment equivalence classes
maxSegId = Seg.Global.getMaxSegId(param);
aLUT = Agglo.buildLUT(maxSegId, aAgglos);
bLUT = Agglo.buildLUT(maxSegId, bAgglos);

aInB = cellfun( ...
    @(ids) setdiff(bLUT(ids), 0), ...
    aAgglos, 'UniformOutput', false);
assert(all(cellfun(@numel, aInB) <= 1));

expIdsA = find(not(cellfun(@isempty, aInB)));
assert(all(arrayfun(@(id) isequal(id, aInB{id}), expIdsA)));

bInA = cellfun( ...
    @(ids) setdiff(aLUT(ids), 0), ...
    bAgglos, 'UniformOutput', false);
assert(all(cellfun(@numel, bInA) <= 1));

expIdsB = find(not(cellfun(@isempty, bInA)));
assert(all(arrayfun(@(id) isequal(id, bInA{id}), expIdsB)));

%% Export example to webKNOSSOS
if ~isempty(nmlDir)
    exportId = find(a.indBigAxons);
    exportId = exportId(nmlId);

    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    skel = skel.addTree( ...
        sprintf('A. Agglomerate %d', exportId), ...
        a.axons(exportId).nodes(:, 1:3), a.axons(exportId).edges);
    skel = skel.addTree( ...
        sprintf('B. Agglomerate %d', exportId), ...
        b.axons(exportId).nodes(:, 1:3), b.axons(exportId).edges);
    
    skel.write(fullfile(nmlDir, sprintf('axon-%d.nml', exportId)));
end
