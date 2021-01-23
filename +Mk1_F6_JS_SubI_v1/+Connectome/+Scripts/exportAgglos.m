% run this after running Connectome.proofreadSnapses.m
% export pre and post agglos as skeletons for proofreading SASD pairs
Util.log('loading graph...')
graph = load(fullfile(rootDir, 'agglomeration/graph.mat'));
graph = graph.graph;

maxSegId = segmentMeta.maxSegId;
points = transpose(segmentMeta.point);

% HACK(amotta): Convert graph table into historical graph table;
graphS = struct;
graphS.edges = graph.edges;

% use idPre and idPost to extract agglos for pre and post
agglosPre = conn.axons(idPre);
agglosPost = conn.dendrites(idPost);

idxOutWrite = 601:1694;

agglosPreOut = agglosPre(idxOutWrite);
agglosPostOut = agglosPost(idxOutWrite);

datasetName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
voxelSize = [11.24, 11.24, 30];

% write nmls
% pre agglos out
for i=1:numel(agglosPreOut)
    outFile = fullfile(rootDir, 'connectome','nmls', 'agglos',['agglo_pre_' num2str(idxOutWrite(1)) '_' num2str(idxOutWrite(end)) '_' num2str(600+i,'%02d')  '.nml']);
    curAgglo = agglosPreOut(i);
    skel = Skeleton.fromAgglo(graphS, points, curAgglo);
    skel = skel.setParams(datasetName, voxelSize, [0, 0, 0]);
    
    description = sprintf('%s (%s)', info.filename, info.git_repos{1}.hash);
    skel = skel.setDescription(description);
    skel.write(outFile);
end

% post agglos out
for i=1:numel(agglosPostOut)
    outFile = fullfile(rootDir, 'connectome','nmls', 'agglos',['agglo_post_' num2str(idxOutWrite(1)) '_' num2str(idxOutWrite(end)) '_' num2str(600+i,'%02d')  '.nml']);
    curAgglo = agglosPreOut(i);
    skel = Skeleton.fromAgglo(graphS, points, curAgglo);
    skel = skel.setParams(datasetName, voxelSize, [0, 0, 0]);

    description = sprintf('%s (%s)', info.filename, info.git_repos{1}.hash);
    skel = skel.setDescription(description);
    skel.write(outFile);
end
