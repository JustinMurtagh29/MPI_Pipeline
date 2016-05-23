function job = interfaceClassificationOnCluster(p, options)
%INTERFACECLASSIFICATIONONCLUSTER Run interface classifier on cluster
% INPUT pFile: Path to segmentation paramter struct.
%       options: (Optional) Struct which can contain the following fields
%           'saveFeatures': Bool specifying whether feature should be
%               saved.
%               (Default: false)
%           'aggloT': see interfaceClassificationForSeg
%           'saveInterfaces': Bool specifying whether interfaces should be
%               saved.
%               (Default: false)
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
if ~isfield(options,'saveFeatures')
    options.saveFeatures = false;
end
if ~isfield(options,'aggloT')
    options.aggloT = [];
end
if ~isfield(options,'saveInterfaces')
    options.saveInterfaces = false;
end


if ~isfield(options,'cubeIndices')
    options.cubeIndices = 1:numel(p.local);
elseif ~isrow(options.cubeIndices)
    options.cubeIndices = options.cubeIndices';
end


inputCell = cell(length(options.cubeIndices),1);
for i = 1:length(options.cubeIndices)
    cubeNo = options.cubeIndices(i);
    inputCell{i} = {p, i, options.aggloT, options.saveFeatures, ...
                    options.saveInterfaces};
end
job = startCPU(@interfaceClassificationForSeg, inputCell, 'interfaceClassification',24);

end
