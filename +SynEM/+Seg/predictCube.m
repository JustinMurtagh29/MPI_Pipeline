function [scores, X, interfaces] = predictCube( p, cubeNo, fm, classifier )
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

pCube = p.local(cubeNo);
fprintf(['[%s] SynEM.Seg.predictCube - Predicting segmentation ' ...
    'cube %s.\n'], datestr(now), pCube.saveFolder);

%load segmentation
seg = loadSegDataGlobal(p.seg, pCube.bboxSmall);

%load svg
m = load(pCube.edgeFile);
edges = m.edges;
m = load(pCube.borderFile);
borders = m.borders;

if size(edges,1) ~= 0 % skip if all zero segmentCube (mirrorPad)
    %calculate interfaces
    interfaces = SynEM.Svg.calculateInterfaces(seg, edges, borders, ...
        fm.areaT, p.raw.voxelSize, fm.subvolsSize);
    if ~isempty(interfaces.surface)
        %load raw
        bboxFM = bsxfun(@plus, pCube.bboxSmall,[-fm.border', fm.border']./2);
        raw = loadRawData(p.raw, bboxFM);
		
        %calculate features
        X = fm.calculate(interfaces, raw);
        
        %classify
        if ischar(classifier)
            m = load(classifier);
            classifier = m.classifier;
        end
        [~,scores] = classifier.predict(X);
    else
        interfaces =  [];
        scores = [];
        X = [];
    end
else
    interfaces =  [];
    scores = [];
    X = [];
end
