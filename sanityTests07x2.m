%% Load both parameter sets
pOld = '/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat';
pNew = '/gaba/u/mberning/results/pipeline/test/allParameter.mat';

%% Check whether classification yields same result
bbox = [1000 1200; 1000 1200; 1000 1200];

classOld = loadClassData(pOld.p.class.root, pOld.p.class.prefix, bbox);
classNew = loadClassData(pNew.p.class.root, pNew.p.class.prefix, bbox);

all(classOld(:) == classNew(:));

%% Check whether segmentation on inner cube yields same result

%% Check correspondences

%% Check globalization

%% Check graph construction

%% Check feature calculation
