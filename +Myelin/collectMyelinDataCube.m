function voxMyelin = collectMyelinDataCube(p,cubeIdx)

bbox = p.local(cubeIdx).bboxSmall;

load(fullfile(p.local(cubeIdx).saveFolder,'segmentMeta.mat'),'segIds');
myelin = loadSegDataGlobal(p.myelin, bbox);
seg = loadSegDataGlobal(p.seg, bbox);


voxMyelin = seg(myelin(:));

voxMyelin = arrayfun(@(x) sum(x==voxMyelin),segIds);

