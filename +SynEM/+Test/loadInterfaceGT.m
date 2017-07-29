function testSet = loadInterfaceGT( dataPath, smallBBox, excOnly, renumber)
%LOADINTERFACEGT Load the synaptic interfaces dense ground truth.
% INPUT dataPath: string
%           Path to SynEM data folder.
%       smallBBox: (Optional) logical
%           Flag indicating whether only interfaces with center of mass at
%           least 160 nm away from the cube borders are considered.
%           (Default: false)
%       excOnly: (Optional) logical
%           Interfaces overlapping with inhibitory synapses are deleted.
%           (Default: false)
%       renumber: (Optional) logical
%           Flag if group should be renumbered starting from 1.
%           (Default: true)
% OUTPUT testSet: struct
%           Struct containing the following fields
%        group: [Nx1] int
%           Synaptic interfaces labeled by the ground truth segmentation.
%           All interfaces with the same id tentatively belong to one
%           synapse.
%        X: [NxM] single
%           The features for the corresponding rows in group.
%        scores: [Nx2] double
%           SynEM predictions for the test dataset.
%        globalEdgeIdx: [Nx1] int
%           The linear global edge indices for the output interfaces.
%        localEdgeIdx: [Nx1] int
%           The linear local edge indices for the output interfaces, i.e.
%           within the local segmentation cube.
%        isInh: [Nx1] logical
%           Logical indices of all interfaces in group that are inhibitory
%        toKeep: [Nx1] logical
%           Logical indices of the final set of interfaces among all
%           interfaces above the area threshold. toKeep indices are
%           affected by the input flags.
%        borders: [Nx1] struct
%           The borders struct for the test set.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('renumber', 'var') || isempty(renumber)
    renumber = true;
end

%get test cube data
m = load(fullfile(dataPath, 'TestSet', 'TestCubeAux.mat'));
testCube = m.testCube;
globalEdgeIdx = testCube.globalEdgeIdx(1):testCube.globalEdgeIdx(2);
globalEdgeIdx = globalEdgeIdx(testCube.areaT_idx);
localEdgeIdx = find(testCube.areaT_idx);

p = SynEM.Test.loadSegParams(dataPath);
pCube = p.local(67);
mTest = load(fullfile(dataPath, 'TestSet', 'SynapseDetectionTestSet.mat'));
toKeep = true(length(mTest.interfaceLabels), 1);

%crop to bbox small
if exist('smallBBox', 'var') && smallBBox
    m = load(pCube.borderFile);
    borders = m.borders(testCube.areaT_idx);
    cen = {borders(:).Centroid}';
    cen = cell2mat(cen);
    cen = uint16(bsxfun(@plus, cen, pCube.bboxSmall(:,1)'-[1, 1, 1]));
    %get inner bbox
    bbox = pCube.bboxSmall;
    border = ceil(160./[11.24; 11.24; 28]);
    bboxInner = bsxfun(@plus,bbox,[border,-border]);
    isInBboxSmall = all(bsxfun(@ge,cen,bboxInner(:,1)'),2) & ...
                all(bsxfun(@le,cen,bboxInner(:,2)'),2);
    toKeep = bsxfun(@and, toKeep, isInBboxSmall);
end

%get inhibitories
isInh = ismember(mTest.interfaceLabels, mTest.inhSynapses);

%load outputs
group = mTest.interfaceLabels;
m = load(fullfile(pCube.saveFolder, 'InterfaceFeatures.mat'));
X = m.X;
m = load(pCube.synapseFile);
scores = reshape(m.scores, [], 2);

if exist('excOnly', 'var') && excOnly
    toKeep = bsxfun(@and,toKeep, ~isInh);
end

%restrict output according to input options
group = group(toKeep);
testSet.X = X(toKeep, :);
testSet.scores = scores(toKeep, :);
testSet.globalEdgeIdx = globalEdgeIdx(toKeep);
testSet.localEdgeIdx = localEdgeIdx(toKeep);
testSet.isInh = isInh(toKeep);
testSet.borders = borders(toKeep);
testSet.toKeep = toKeep;

%get rid of empty groups
if renumber
    group(group > 0) = Util.renumber(group(group>0), 'stable'); 
end
testSet.group = group;
end

