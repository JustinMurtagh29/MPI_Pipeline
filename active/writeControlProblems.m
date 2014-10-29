function writeControlProblems( pT )
% Write control problems based on skeletons

display('Loading data:');
tic;
% load data
segNew = loadSegData(pT.seg.root, pT.seg.prefix, pT.bboxBig);
load(pT.edgeFile);
load(pT.borderFile);
train = load([pT.gtSkel '01.mat']);
toc

display('Processing segments according to ground truth:');
tic;
% Align training edges with border structure
isConnected = zeros(size(edges,1),1);
for i=1:size(edges,1)
    % Compare only with skeleton training data
    edgeFound = all(bsxfun(@eq, train.edges, edges(i,:)),2);
    if any(edgeFound)
        allLabels = train.labels(edgeFound);
        isConnected(i) = unique(allLabels);
    end
end
% Replace classification with ground truth
edgesOld = edges;
pOld = (isConnected + 1)./2;
borderOld = border;
% Remove edges with probability below lower cutoff from query (not from possible neighbours)
rejected = pOld < 0.75;
edgesNew = edgesOld(~rejected,:);
pNew = pOld(~rejected);
borderNew = borderOld(~rejected);
% Remove all borders outside inner bbox from query (same)
lowerBounds = -pT.border(:,1)';
upperBounds = pT.bboxSmall(:,2)' - pT.bboxSmall(:,1)' - pT.border(:,1)';
lowerBounds = lowerBounds([2 1 3]); % this switch is annoying, any nicer solution, continues in writeKnowledgeDB?
upperBounds = upperBounds([2 1 3]);
outsideBbox = any(bsxfun(@lt,cat(1,borderNew(:).Centroid),lowerBounds),2) | any(bsxfun(@gt,cat(1,borderNew(:).Centroid),upperBounds),2);
edgesNew = edgesNew(~outsideBbox,:);
pNew = pNew(~outsideBbox);
borderNew = borderNew(~outsideBbox);
toc

% Do the same as in fromGraphToDB otherwise
display('Mitochondria detetcion:');
tic;
% Mitochondria detektieren
raw = loadRawData(pT.raw.root, pT.raw.prefix, pT.bboxBig, 0); 
mito = mitoDetection(raw, segNew);
mito = uint8(mito);
toc

display('Filling watershed lines (dilation):');
tic;
% Segmentation auswachsen
segTemp = imdilate(segNew, ones(3,3,3));
borders = segNew == 0;
segNew(borders) = segTemp(borders);
toc

display('writing to knowledge DB:');
tic;
% transfer everything to knowledgeDB
pT.controlFlag = 1;
save(['/zdata/manuel/sync/problemInspector/' 'kDbBefore' strrep(pT.start, '/', '-') '.mat'], '-v7.3');
missions = writeKnowledgeDB(pT, segNew, mito, edgesNew, pNew, borderNew, edgesOld, pOld, borderOld);
toc

end
