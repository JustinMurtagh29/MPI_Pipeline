% script to evaluate the axon agglomeration using the new GT axons
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();

ovT = 2; % minimal overlap in segments

%% load axon agglo

p = Gaba.getSegParameters('ex145_ROI2017');
p.agglo.axonAggloFile = fullfile(p.agglo.saveFolder, ...
    'axons_19_a_partiallySplit_v2.mat');
axons = L4.Axons.getLargeAxons(p, false, false);


%% get the overlaps between gt skeletons and agglos


[skels, segIds] = connectEM.eval.getNewAxonGT();
for i = 1:10
    skels{i}.verbose = false;
end

m = load(p.agglo.axonAggloFile);
ov = connectEM.eval.getNewAxonGTAggloOverlap(segIds, axons, ovT);

[~, axFile] = fileparts(p.agglo.axonAggloFile);
outFile = fullfile(p.agglo.saveFolder, 'eval', ...
    sprintf('axon_gt_eval_%s.mat', axFile));

Util.ssave(outFile, info, ov, skels);


%% everything to nml

if iscell(axons) % if agglos are used then create MST superagglo
    m = load(p.svg.segmentMetaFile, 'point');
    point = point';
    axons = SuperAgglo.fromAgglo(axons, point, 'mst', 'voxelSize', ...
        [11.24, 11.24, 28]);
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
    
    [skel, gid] = skel.addGroup(sprintf('Axon_%02d', i));
    skel = skel.addTreesToGroup(c:c+numT, gid);
    c = c + numT;
end

skel = skel.setDescription(sprintf(['Overlap (ovT=%d) of axon agglo %s ' ...
    'with axon_gt_new.'], ovT, axFile));
skel = Skeleton.appendRunInfoToDescription(skel, info);
skel.write(sprintf('AxonGT_%s_overlaps.nml', axFile));
