function [data,metadata] = calculateInterfacesFromSeg( p, cubeNo, featureMap, idConversion, cubeInfo )
%CALCULATEINTERFACESFROMSEG Calculate interfaces using Manuels segmentation
%and interface calculations.
% INPUT p: Segmentation parameter file.
%       cubeNo: Number of cube to process as linear index (in p.local).
%       featureMap: The feature map used for the feature calculation. Must
%           be a struct containing the fields
%            - areaThreshold: area threshold for interfaces
%            - rinclude: vector with size of subsegments
%            - border: Border around the inner bounding box for later
%                      feature calculation.
%       idConversion: (Optional) Convert neighborIDs to local ids if
%            required. Local IDs are necessary here. Use 'GlobalToLocal' to
%            convert to local indices or 'none'.
%            (Default: 'none')
%       cubeInfo: For synapse prediction only specify a single cube from
%                 p.local, set the first two arguments to [] and only
%                 consider the first output argument.
% OUTPUT data: Struct containing the fields
%           - raw: Raw data array with coordinates metadata.bboxBig
%                  containing a the metadata.border around the region where
%                  interfaces are calculated from the dataset
%                  metadata.experiments
%           - segments: Segments array in which the interfaces are
%                       calculated with coordinates metadata.bboxSmall from
%                       the segmentation metadata.segmentation in cube
%                       metadata.cubeNo
%           - interfaceSurfaceList: Linear coordinates of interface voxels
%                       segments with at least metadata.areaThreshold
%                       voxels
%           - subsegmentsList: Subsegments calculated in segments for
%                       different metadata.rinclude.
%           - neighborIDs: segment IDs for corresponding interfaces in
%                       coordinates specified in idConversion (none = local
%                       ids)
%           - centroidList: Center of mass/centroids of corresponding
%                       interfaces.
%        metadata: Struct containing the information mentioned in data.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%load interfaces
if ~exist('idConversion','var') || isempty(idConversion)
    idConversion = 'none';
end
if ~exist('cubeInfo','var') || isempty(cubeInfo)
    cubeInfo = p.local(cubeNo);
end
[ interfaceSurfaceList, centroidList, neighborIDs ] = loadInterfacesFromSeg( cubeInfo, featureMap.areaThreshold, idConversion );

%load segmentation
m = matfile(cubeInfo.segFile);
segments = m.seg;
from = cubeInfo.bboxSmall(:,1) - cubeInfo.bboxBig(:,1) + 1;
to = from + cubeInfo.bboxSmall(:,2) - cubeInfo.bboxSmall(:,1);
segments = segments(from(1):to(1),from(2):to(2),from(3):to(3));

subsegmentsList = calculateSubsegments(interfaceSurfaceList,neighborIDs,segments,featureMap.rinclude, p.raw.voxelSize);

%load raw and make output
bboxBig = cubeInfo.bboxSmall + [-featureMap.border,featureMap.border];
raw = readKnossosRoi(p.raw.root,p.raw.prefix,bboxBig,'uint8');
if all(~raw)
    error('The raw file read from %s is empty.',knossos_conf.parfolder);
end

%generate output
data.raw = raw;
data.segments = segments;
data.interfaceSurfaceList = interfaceSurfaceList;
data.subsegmentsList = subsegmentsList;
data.neighborIDs = neighborIDs;
data.centroidList = centroidList;

if nargout == 2
    metadata.experiment = p.raw.prefix;
    metadata.segmentation = p.start;
    metadata.bboxSmall = cubeInfo.bboxSmall;
    metadata.bboxBig = bboxBig;
    metadata.areaThreshold = featureMap.areaThreshold;
    metadata.rinclude = featureMap.rinclude;
    metadata.idConversion = idConversion;
    metadata.cubeNo = cubeNo;
    metadata.border = featureMap.border;
end

end

