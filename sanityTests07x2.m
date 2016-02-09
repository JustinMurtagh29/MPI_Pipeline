%% Load both parameter sets
pOld = load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
pNew = load('/gaba/u/mberning/results/pipeline/test/allParameter.mat');

%% Check whether classification yields same result
bbox = [1200 1400; 1200 1400; 1200 1400];

classOld = loadClassData(pOld.p.class.root, pOld.p.class.prefix, bbox);
classNew = loadClassData(pNew.p.class.root, pNew.p.class.prefix, bbox);

sum(abs(classOld(:) - classNew(:)) < 1e-5)
numel(classOld)

% Conclusion, close enough see changes matlab R2011b -> R2015b for explanantion of 2 order of magnitude float point eps 1e-7

%% Check whether segmentation on inner cube yields same result

segOld = load(pOld.p.local(2,2,2).segFile);
segNew = load(pNew.p.local(2,2,2).tempSegFile);

sum(segOld.seg(:) == segNew.seg(:))
numel(segOld.seg)

% Same conclusion

%% Check correspondences between to inner cubese

% Due to symbolic linking in ealier pipelines (change in FOV): 

cOld = load('/gaba/u/mberning/results/pipeline/20141007T094904/correspondences/040403050403.mat');
cNew = load([pNew.p.correspondence.saveFolder '030403040403.mat']);

cNew.uniqueCorrespondences == bsxfun(@plus, single(cOld.result.correspondences), [0 1])

% Conclusion: One new segment in 'upper' cube, thus shift of one in correspondences, good otherwise

%% Check cutting out of inner bounding box

segOld = load(pOld.p.local(2,2,2).segFile);
segOld.seg = segOld.seg(257:end-256,257:end-256,129:end-128);
segNew = load(pNew.p.local(2,2,2).segFile);

sum(segOld.seg(:) == segNew.seg(:))
numel(segOld.seg)

% Still the same here
unique(segNew.seg(:))
unique(segOld.seg(:))

%% Check globalization

% Check whether segmentation in segGlobal is consitent (other global IDs now, continous)
segOld = load([pOld.p.local(2,2,2).saveFolder 'segGlobal.mat']);
segOld.seg = segOld.seg(257:end-256,257:end-256,129:end-128);
segNew = load([pNew.p.local(2,2,2).saveFolder 'segGlobal.mat']);

% Check whether borders are at same position (- noise, see bigFwdPass)
sum((segOld.seg(:) > 0) == (segNew.seg(:) > 0))
numel(segOld.seg)

% NOT the same here
unique(segNew.seg(:))
unique(segOld.seg(:))

% Check whether correspondences (globalized) are consitent
cOld = load('/gaba/u/mberning/results/pipeline/20141007T094904/correspondences/040403050403global.mat');
cNew = load([pNew.p.correspondence.saveFolder '030403040403global.mat']);

cNew.uniqueCorrespondences == bsxfun(@plus, single(cOld.result.correspondences), [0 1])

%% Check graph construction



%% Check feature calculation



