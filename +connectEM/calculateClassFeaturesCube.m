function calculateClassFeaturesCube(p, cubeNo, fm)
% Calculate SynEM features for SegEM classification output on a local cube

% load segmentation
seg = loadSegDataGlobal(p.seg, p.local(cubeNo).bboxSmall);

%load svg
load(p.local(cubeNo).edgeFile, 'edges');
load(p.local(cubeNo).borderFile, 'borders');

%calculate interfaces
interfaces = SynEM.Svg.calculateInterfaces(seg, edges, borders, ...
    fm.areaT, p.raw.voxelSize, fm.subvolsSize);

%load raw
bboxFM = bsxfun(@plus, p.local(cubeNo).bboxSmall,[-fm.border', fm.border']./2);
class = loadClassData(p.class, bboxFM);

%calculate features
features = fm.calculate(interfaces, class);

% Save features
if strcmp(fm.mode, 'direction')
    features = features(1:end/2,:); %only save first direction
end
outputFile = [p.local(cubeNo).saveFolder 'InterfaceClassFeatures.mat'];
Util.save(outputFile, features);

end
