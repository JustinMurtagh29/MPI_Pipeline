function fromGraphToDB( p )

display('Loading data:');
tic;
% load data
seg = loadSegData(p.seg.root, p.seg.prefix, p.bboxBig);
load(p.edgeFile);
load(p.borderFile);
load(p.hyper);
train = loadTrainingData(p);
test = loadTestData(p);
toc

display('GP classification:');
tic;
% gpml toolbox usage
run('/zdata/manuel/code/active/gpml/startup.m');
% Make predictions
[a b c d lp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, train.X, train.y, test.X, ones(size(test.X,1), 1));
toc

display('Joining segments according to classification results:');
tic;
% join segments according to probability cutoffs given in paramBG
prob = exp(lp);
lowerBounds = -p.border(:,1)';
upperBounds = p.bboxSmall(:,2)' - p.bboxSmall(:,1)' - p.border(:,1)';
[segNew, edgesNew, pNew, borderNew, edgesOld, pOld, borderOld] = joinSegments(seg, edges, prob, border, p.lowerCut, p.upperCut, lowerBounds, upperBounds);
toc

display('Mitochondria detetcion:');
tic;
% Mitochondria detektieren
raw = loadRawData(p.raw.root, p.raw.prefix, p.bboxBig); 
raw = single(raw);
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

display('Saving data and writing to knowledge DB:');
tic;
% transfer everything to knowledgeDB, resort first
[pNew, idx] = sort(pNew, 'descend');
edgesNew = edgesNew(idx,:);
borderNew = borderNew(idx);
save(['/zdata/manuel/sync/problemInspector/' 'kDbBefore' strrep(p.start, '/', '-') '.mat'], '-v7.3');
%missions = writeKnowledgeDB(p, segNew, mito, edgesNew, pNew, borderNew, edgesOld, pOld, borderOld);
toc

end

