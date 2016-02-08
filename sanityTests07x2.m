%% Load both parameter sets
pOld = load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
pNew = load('/gaba/u/mberning/results/pipeline/test/allParameter.mat');

%% Check whether classification yields same result
bbox = [1200 1400; 1200 1400; 1200 1400];

classOld = loadClassData(pOld.p.class.root, pOld.p.class.prefix, bbox);
classNew = loadClassData(pNew.p.class.root, pNew.p.class.prefix, bbox);

sum(abs(classOld(:) - classNew(:)) < 1e-5)
numel(classOld)

%% Check whether segmentation on inner cube yields same result

%% Check correspondences between to inner cubese

%% Check globalization

%% Check graph construction

%% Check feature calculation
