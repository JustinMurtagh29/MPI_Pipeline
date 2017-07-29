function [ X, interfaces ] = cubeFeatures( p, cubeNo, fm )
%CUBEFEATURES Calculate the interfaces features for a local segmentation
%cube.
% INPUT p: struct
%           SegEM segmentation parameter struct.
%       cubeNo: int
%           Linear index of the cube in p.local
%       fm: SynEM.FeatureMap
%           The feature map for prediction.
%       classifier: object or string
%           Classifier object (e.g. SynEM.Classifier) or path to classifier
%           object mat-file containing a 'classifier' variable.
% OUTPUT X: [NxM] single
%           The feature map used for prediction. Rows correspond to
%           interfaces, columns to features.
%        interfaces: struct
%           Voxel data of each interface. Fields are
%           surface: [Nx1] cell array where each cell contains the linear
%               indices of voxels w.r.t. seg of an interface surface as an
%               [Nx1] int array.
%           subseg: [1xN] cell array where N = length(rinclude). Each cell
%               contains a [Mx2] cell array where
%               M = length(interfaceSurface) which contains the linear
%               indices of the first and second subsegment w.r.t. seg.
%
% NOTE Only borders with at least fm.areaT voxels are considered and thus
%      size(scores,1) is length number of borders > fm.areaT.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

pCube = p.local(cubeNo);
fprintf(['[%s] SynEM.Seg.predictCube - Predicting segmentation ' ...
    'cube %s.\n'], datestr(now), pCube.saveFolder);

%load segmentation
seg = SynEM.Aux.readKnossosRoi(p.seg.root, p.seg.prefix, ...
    pCube.bboxSmall, 'uint32', '', 'raw');

%load svg
m = load(pCube.edgeFile);
edges = m.edges;
m = load(pCube.borderFile);
borders = m.borders;

%calculate interfaces
interfaces = SynEM.Svg.calculateInterfaces(seg, edges, borders, ...
    fm.areaT, p.raw.voxelSize, fm.subvolsSize);

%load raw
bboxFM = bsxfun(@plus, pCube.bboxSmall,[-fm.border', fm.border']./2);
raw = SynEM.Aux.readKnossosRoi(p.raw.root, p.raw.prefix, bboxFM);

%calculate features
X = fm.calculate(interfaces, raw);

end

