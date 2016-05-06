function job = calculateFeaturesForSegOnCluster( pFile, featureMap, options )
%CALCULATEFEATURESFORSEGONCLUSTER Run interface feature calculation on
% cluster.
% INPUT pFile: Path to segmentation paramter struct.
%       featureMap: An interface classifier feature map.
%       options: (Optional) Struct which can contain the following fields
%           'saveInterfaces': Bool specifying whether interfaces should be
%               saved.
%               (Default: false);
%           'aggloT': see calculateFeaturesForSeg
%           'cubeIndices': [Nx1] array of integer specifying the cubes in
%               p.local for which the features are calculated.
%               (Default: all cubes)
% OUTPUT job: Matlab job object.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse options
if ~exist('options','var') || isempty(options)
    options = struct();
end
if ~isfield(options,'saveInterfaces')
    options.saveInterfaces = false;
end
if ~isfield(options,'aggloT')
    options.aggloT = [];
end

%load parameter file
m = load(pFile);
p = m.p;

if ~isfield(options,'cubeIndices')
    options.cubeIndices = 1:numel(p.local);
elseif ~isrow(options.cubeIndices)
    options.cubeIndices = options.cubeIndices';
end

cluster = getCluster('cpu');
inputCell = cell(length(options.cubeIndices),1);
for i = 1:length(options.cubeIndices)
    cubeNo = options.cubeIndices(i);
    inputCell{i} = {pFile, cubeNo, featureMap, options.aggloT, ...
                    options.saveInterfaces};
end
job = startJob(cluster,@calculateFeaturesForSeg, inputCell, 0);


end

