function calculateFeaturesForSeg( pFile, cubeNo, featureMap, aggloT, saveInterfaces )
%CALCULATEFEATURESFORSEG Calculation of interfaces and features.
% INPUT pFile: Path to segmentation paramter struct.
%       cubeNo: Integer specifying the linear index of the local segmentation
%           cube in p.local.
%       featureMap: See makeFeatureMap.
%       aggloT: (Optional) Double between 0 and 1 specifying the lower
%           threshold on the GP probability above which segments will get
%           merged.
%           (Default: 1 - no prior merging).
%       saveInterfaces: (Optional) Logical specifying whether to save the
%           interfaces.
%           Interfaces are saved to
%           [p.local(i).saveFolder 'interfaces.mat'].
%           (Default: false)
%
% see also interfaceClassification
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~exist('saveInterfaces','var') || isempty(saveInterfaces)
    saveInterfaces = false;
end

%load parameter file
m = load(pFile);
p = m.p;
pCube = p.local(cubeNo);

%load segmentation
seg = Seg.Local.getSegSmall(pCube, true);

%load edges and borders
m = load(pCube.edgeFile);
edges = m.edges;
m = load(pCube.borderFile);
borders = m.borders;

%load raw
bboxRaw = pCube.bboxSmall + [-featureMap.border,featureMap.border];
raw = readKnossosRoi(p.raw.root,p.raw.prefix,bboxRaw,'uint8');
if all(~raw)
    error('The raw file read from %s is empty.',p.raw.root);
end

%prior agglomeration
if exist('aggloT','var') && ~isempty(aggloT)
    m = load(pCube.probFile);
    prob = m.prob;
    m = load(pCube.segmentFile);
    segments = m.segments;
    eClasses = Graph.findConnectedComponents(edges(prob > aggloT,:));
    [ seg, edges, borders ] = Seg.Local.applySegEquiv( eClasses, seg, edges, borders, segments );
end

%interface and feature calculation
interfaces = calculateInterfaces(seg, edges, borders, featureMap.areaThreshold, featureMap.voxelSize, featureMap.rinclude);
X = calculateFeaturesOld(raw, interfaces, featureMap);

%save features
m = matfile([p.local(cubeNo).saveFolder 'interfaceWeights.mat'], 'Writable', true);
m.X = X;

%save additional information in case of agglomeration
if exist('aggloT','var') && ~isempty(aggloT)
    m.aggloT = aggloT;
    tmp = [borders(:).Area] > 150;
    m.edges = edges(tmp,:);
    borders = borders(tmp);
    m.edgeCentroids = reshape([borders(:).Centroid], 3, length(borders))';
end
    
if saveInterfaces
    m = matfile([p.local(cubeNo).saveFolder 'interfaces.mat'], 'Writable', true);
    m.interfaces = interfaces;
end

end

