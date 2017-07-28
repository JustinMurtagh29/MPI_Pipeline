function [scores, X, interfaces] = predictCube( p, cubeNo, fm, ...
    classifier, normalizeRaw )
%PREDICTCUBE Prediction for one segmentation cube.
% INPUT p: struct
%           SegEM segmentation parameter struct.
%       cubeNo: int
%           Linear index of the cube in p.local
%       fm: SynEM.FeatureMap
%           The feature map for prediction.
%       classifier: object or string
%           Classifier object (e.g. SynEM.Classifier) or path to classifier
%           object mat-file containing a 'classifier' variable.
%       normalizeRaw: (Optional) logical
%           Flag to indicate that the raw data is normalized to 122 mean
%           and 22 standard deviation by using p.norm.func to normalize to
%           mean 0 and 1 std first.
%           (Default: no raw data normalization).
% OUTPUT scores: [Nx1] double
%           Prediction score for each interface in the local segmentation
%           cube.
%        X: [NxM] single
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
%      If the 'direction' mode is used than each interface produces two
%      scores (one for each direction) which is arranged such that the
%      first half of the scores is one direction for each interface and the
%      second half of the scores is the other direction (i.e. typically
%      reshape(scores,[],2) gives the two scores for one interface in the
%      same row).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('normalizeRaw','var') || isempty(normalizeRaw)
    normalizeRaw = false;
end

pCube = p.local(cubeNo);
Util.log('Predicting segmentation cube %s.', pCube.saveFolder);

%load segmentation
seg = SynEM.Aux.readKnossosRoi(p.seg.root, p.seg.prefix, ...
    pCube.bboxSmall, 'uint32', '', 'raw');

% load svg
m = load(pCube.edgeFile);
edges = m.edges;
m = load(pCube.borderFile);
borders = m.borders;

% stop here if there are no borders above the area threshold
if ~any([borders.Area] > fm.areaT)
    scores = zeros(0, 1);
    X = [];
    interfaces.surface = [];
    interfaces.subseg = [];
    return
end

% calculate interfaces
interfaces = SynEM.Svg.calculateInterfaces(seg, edges, borders, ...
    fm.areaT, p.raw.voxelSize, fm.subvolsSize);

% load raw
bboxFM = bsxfun(@plus, pCube.bboxSmall, [-fm.border', fm.border']./2);
raw = SynEM.Aux.readKnossosRoi(p.raw.root, p.raw.prefix, bboxFM);

% raw data normalization
if normalizeRaw
    Util.log('Normalizing raw data to ex145 statistics.');
    raw = double(raw);
    raw = p.norm.func(raw);
    raw = single(raw.*22 + 122);
end

% calculate features
X = fm.calculate(interfaces, raw);

% classify
if ischar(classifier);
    m = load(classifier);
    classifier = m.classifier;
end
[~,scores] = classifier.predict(X);

end
