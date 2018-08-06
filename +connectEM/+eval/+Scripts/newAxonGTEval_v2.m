% script to evaluate the axon agglomeration using the new GT axons
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


ovT = 2; % minimal overlap in segments
nhood = 0;

info = Util.runInfo();


%% load axon agglo

p = Gaba.getSegParameters('ex145_ROI2017');
p.agglo.axonAggloFile = fullfile(p.agglo.saveFolder, ...
    'axons_19_a_partiallySplit_v2.mat');
[axons, ~, ~, axFile] = L4.Axons.getLargeAxons(p, true, true);


%% get the overlaps between gt skeletons and agglos


[skels, segIds] = connectEM.eval.getNewAxonGT([], nhood);
for i = 1:10
    skels{i}.verbose = false;
end

ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, axons, ovT);
[stats, debug] = connectEM.eval.gtReconstructionStats(skels, segIds, axons, ov);

[~, axFile] = fileparts(axFile);
outFile = fullfile(p.agglo.saveFolder, 'eval', ...
    sprintf('axon_gt_eval_%s.mat', axFile));

Util.ssave(outFile, info, ov, skels);


%% everything to nml

if iscell(axons) % if agglos are used then create MST superagglo
    m = load(p.svg.segmentMetaFile, 'point');
    point = m.point';
    axonsBkp = axons;
    idx = unique(cell2mat(cellfun(@(x)x(:,1), ov, 'uni', 0)));
    tmp = SuperAgglo.fromAgglo(axonsBkp(idx), point, 'mst', ...
        'voxelSize', ...
        [11.24, 11.24, 28]);
    clear axons
    axons(idx) = tmp;
end

ovSkels = connectEM.eval.newAxonGTOverlapsToSkel(skels, axons, ov);

% write to tracings
skel = L4.Util.getSkel();
c = 1;
for i = 1:10
    toMergeSkel = ovSkels{i};
    numT = toMergeSkel.numTrees();
    
    skel = skel.mergeSkels(toMergeSkel);
    skel.colors{c} = [0, 1, 0, 1];
    skel.colors(c+1:end) = {[1, 0, 0, 1]};
    skel.names{c} = [skel.names{c}, ...
        sprintf(' (path length: %.2f um; recall: %.3f)', ...
        stats.pathLength(i) ./ 1000, stats.recall(i))];
    
    [skel, gid] = skel.addGroup(sprintf('Axon_%02d', i));
    skel = skel.addTreesToGroup(c:c+numT, gid);
    c = c + numT;
end

skel = skel.setDescription(sprintf(['Overlap (ovT=%d) of axon agglo %s ' ...
    'with axon_gt_new.'], ovT, axFile));
skel = Skeleton.appendRunInfoToDescription(skel, info);
% skel.write(sprintf('AxonGT_%s_overlaps.nml', axFile));
