
% Bounding box of _demo in new2
bbox(:,1) = [1597; 3837; 1412];
bbox(:,2) = bbox(:,1) + ceil(20000./[11.24; 11.24; 28]);

root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2017/segmentation/1/';
prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';
data = readKnossosRoi(root, prefix, bbox, 'uint32');
writeKnossosRoi('/tmpscratch/mberning/segForDemo/1/', '2012-09-28_ex145_07x2_demo_mag1', [1 1 1], data, 'uint32');

% Create hierachy for zooming out etc.
bbox = [1,1,1;(ceil(ceil([1; 1; 1] + 20000./[11.24; 11.24; 28])./1024).*1024)']';
createResolutionPyramid('/tmpscratch/mberning/segForDemo/1/', '2012-09-28_ex145_07x2_demo_mag1', bbox, '/tmpscratch/mberning/segForDemo/', true);

bbox(:,1) = p.local(4,9,6).bboxSmall(:,1);
bbox(:,2) = p.local(6,11,7).bboxSmall(:,2);
raw = loadRawData(p.raw, bbox);
class = loadClassData(p.class, bbox);
seg = loadSegDataGlobal(p.seg, bbox);
segmentMeta = rmfield(segmentMeta, {'box' 'centroid'});
graph = rmfield(graph, {'neighbours' 'neighBorderIdx' 'neighProb' 'borderIdx'});

