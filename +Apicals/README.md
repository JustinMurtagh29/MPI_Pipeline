See my *lablog post(s)* for more information.
---

And you also need to have the Agglo Folder from Alessandros repo in your path,
at least the function for calculating path length because the AIS script is dependent of it:
	`Agglo.calcPathLengths(p, aggloComponents);`


---------------------------------------------------------------------------
Use this for reproduction if you are not using the scripts:
```
dendAggloStatePath = '/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_03.mat';
spineHeadsStatePath = '/+Apicals/importantFiles/20170829_spineheadsAttached.mat';
dendLensPath = '/+Apicals/importantFiles/20170829_DendriteAggloLengths.mat';
parameterPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat';
metaPath = '/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat';
local = 'yourmointpoint' %local = '/mnt/gaba';
somaStatePath = '/gaba/u/rhesse/forBenedikt/somasNoClosedHoles.mat';
spineHeadCountsPath = '/+Apicals/importantFiles/20170829_spineheadCounts(OnlyBig).mat';
writePath = '/somePath/20170830_funTest_dthr4000.nml';
```
