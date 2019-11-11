% find out why the agglomeration stopped while growing
% load p, graph

outputFolder = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/20191107T162813_agglomeration/dendrites_border_50seg_100prob_0.99size_1000/'

segId = 16879787

edges= graph.edges;
prob = graph.prob;
borderIdx = graph.borderIdx;
com = segmentMeta.point;

allEdges= ismember(edges,segId);
idxFound = find(allEdges(:,1) | allEdges(:,2));
outEdges = edges(idxFound,:);
outProb = prob(idxFound,:);
outBorderIdx = borderIdx(idxFound);

skel = skeleton();

skel.parameters.experiment.name=p.experimentName;
skel.parameters.scale.x = num2str(p.raw.voxelSize(1));
skel.parameters.scale.y = num2str(p.raw.voxelSize(2));
skel.parameters.scale.z = num2str(p.raw.voxelSize(3));
skel.parameters.offset.x = '0';
skel.parameters.offset.y = '0';
skel.parameters.offset.z = '0';

theseCoM = com(outEdges,:);
theseEdgesSegId = outEdges;

% remove duplicate edges
[theseEdgesSegIds, idxU] = unique(theseEdgesSegId, 'rows');
theseProb = outProb(idxU);
theseEdgesNodes = changem(double(theseEdgesSegIds), 1:size(theseCoM,1), outEdges(:));
theseBorderIdx = outBorderIdx(idxU);

for i=1:size(theseEdgesNodes,1)
    curEdges = theseEdgesNodes(i,:);
    curNodes = theseCoM(curEdges,:);
    curEdgesSegIds = theseEdgesSegIds(i,:);
 
    curProb = theseProb(i);
    curSegId = setdiff(curEdgesSegIds(:), segId);
    dendProb = segmentMeta.dendriteProb(curSegId);
    axonProb = segmentMeta.axonProb(curSegId);
    segSize = segmentMeta.voxelCount(curSegId);
    if isnan(theseBorderIdx(i))
        curBorderSize = NaN;
    else
        curBorderSize = borderMeta.borderSize(theseBorderIdx(i));
    end

    skel = skel.addTree(['tree_conn:' num2str(curProb,'%.3f') '_dend:' num2str(dendProb,'%.3f') '_axon:' num2str(axonProb,'%.3f') '_vx:' num2str(segSize)  '_border:' num2str(curBorderSize)], curNodes);
end

skel.write(fullfile(outputFolder, ['debugSegId' num2str(segId) '.nml']))

