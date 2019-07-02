% find out if spines are being attached in the graph partition

%segId = 7222660; % segment ID of spinehead
segId = 7011483
% find which agglo contains spine head id
maxSegId = Seg.Global.getMaxSegId(p);
idxFound = Agglo.findParentAggloForThisSeg(segId, agglos, maxSegId);

allEdges= ismember(edges,segId);
outEdges = edges(find(allEdges(:,1) | allEdges(:,2)),:);
outProb = prob(find(allEdges(:,1) | allEdges(:,2)),:);

skel = skeleton();

skel.parameters.experiment.name=p.experimentName;
skel.parameters.scale.x = num2str(p.raw.voxelSize(1));
skel.parameters.scale.y = num2str(p.raw.voxelSize(2));
skel.parameters.scale.z = num2str(p.raw.voxelSize(3));
skel.parameters.offset.x = '0';
skel.parameters.offset.y = '0';
skel.parameters.offset.z = '0';

com = segmentMeta.point;
theseCoM = com(outEdges,:);
theseEdgesSegId = outEdges;
% remove duplicate edges
[theseEdgesSegId, idxU] = unique(theseEdgesSegId, 'rows');
theseProb = outProb(idxU);
theseEdgesNodes = changem(double(theseEdgesSegId), 1:size(theseCoM,1), outEdges(:));

for i=1:size(theseEdgesNodes,1)
curEdges = theseEdgesNodes(i,:);
curNodes = theseCoM(curEdges,:);
skel = skel.addTree(['tree_' num2str(curEdges(1)) '-' num2str(curEdges(2)) '_' num2str(theseProb(i),'%03f') ], curNodes)

end

skel.write(fullfile(outputFolder, ['debugSpine' num2str(segId) '.nml']))

