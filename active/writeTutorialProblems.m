function writeTutorialProblems( pT )

display('Loading data:');
tic;
% load data
seg = loadSegData(pT.seg.root, pT.seg.prefix, pT.bboxBig);
load(pT.edgeFile);
load(pT.borderFile);
toc

display('Processing segments according to ground truth:');
tic;
% join segments according to probability cutoffs given in paramBG
segNew = seg;
skel.file = pT.trainFile;
skel.bbox = pT.bboxSmall;
[labelIdx, labels] = extractGroundTruthFromNml(seg, edges, skel);
edgesOld = edges(labelIdx,:);
pOld = (labels' + 1)/2;
borderOld = border(labelIdx);
% Remove edges with probability below lower cutoff from query (not from possible neighbours)
rejected = pOld < pT.lowerCut;
edgesNew = edgesOld(~rejected,:);
pNew = pOld(~rejected);
borderNew = borderOld(~rejected);
% Remove all borders outside inner bbox from query (same)
lowerBounds = -pT.border(:,1)';
upperBounds = pT.bboxSmall(:,2)' - pT.bboxSmall(:,1)' - pT.border(:,1)';
outsideBbox = any(bsxfun(@lt,cat(1,borderNew(:).Centroid),lowerBounds),2) | any(bsxfun(@gt,cat(1,borderNew(:).Centroid),upperBounds),2);
edgesNew = edgesNew(~outsideBbox,:);
pNew = pNew(~outsideBbox);
borderNew = borderNew(~outsideBbox);
% Use a random subset of 30 control problems
randIdx = randperm(length(pNew), 30);
edgesNew = edgesNew(randIdx,:);
pNew = pNew(randIdx);
borderNew = pNew(randIdx);
toc

display('Mitochondria detetcion:');
tic;
% Mitochondria detektieren
raw = loadRawData(pT.raw.root, pT.raw.prefix, pT.bboxBig, 0); 
mito = mitoDetection(raw, seg);
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
p.tutorialFlag = 1;
missions = writeKnowledgeDB(pT, segNew, mito, edgesNew, pNew, borderNew, edgesOld, pOld, borderOld);
toc

end
