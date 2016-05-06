function [ data_sa, metadata, synapsePrediction, data ] = convertSegForSynapseAnnotator( p, cubeNo, featureMap, idConversion  )
%CONVERTSEGFORSYNAPSEANNOTATOR Calculate interface data for the synapse
% annotator gui.
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
% OUTPUT data_sa: A data struct which can be loaded in the synapse
%                 annotator (subsegments list has one rinclude only to
%                 decrease file size). Note that this struct is not
%                 suitable for feature calculation. Use the data output
%                 instead.
%        metadata: Struct containing information about the data stored in
%                  data.
%        synapsePrediction: Labels for the interfaces if a prediction was
%                         found in the segmentation folder.
%        data: A data struct for the feature calculation(i.e. with the
%              correct border and all subsegment indices).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cubeInfo = p.local(cubeNo);
if ~exist('idConversion','var') || isempty(idConversion)
    idConversion = 'none';
end

%calculate standard data and metadata
[data,metadata] = calculateInterfacesFromSeg( p, cubeNo, featureMap, idConversion );

%get raw and segments for synapse annotator
m = matfile(cubeInfo.segFile);
segments = m.seg;
segments = segments(257 - 100:768 + 100,257 - 100:768 + 100,129 - 40:384 + 40);
raw = readKnossosRoi(p.raw.root,p.raw.prefix,cubeInfo.bboxBig,'uint8');
if all(~raw)
    error('The raw file read from %s is empty.',p.raw.root);
end
raw = raw(257 - 100:768 + 100,257 - 100:768 + 100,129 - 40:384 + 40);

%transform indices in surface list to larger bounding box
interfaceSurfaceList = cell(length(data.interfaceSurfaceList),1);
subsegmentsList = cell(1);
for i = 1:length(data.interfaceSurfaceList)
    indices = data.interfaceSurfaceList{i};
    [x,y,z] = ind2sub([512 512 256],indices);
    x = x + 100;
    y = y + 100;
    z = z + 40;
    indices = sub2ind(size(raw),x,y,z);
    interfaceSurfaceList{i} = uint32(indices);
    for j = 1:2
        indices = data.subsegmentsList{1}{i,j};
        [x,y,z] = ind2sub([512 512 256],indices);
        x = x + 100;
        y = y + 100;
        z = z + 40;
        indices = sub2ind(size(raw),x,y,z);
        subsegmentsList{1}{i,j} = uint32(indices);
    end
end

%get synapse prediction if available
try
    m = matfile(cubeInfo.synapseFile);
    scores = m.scores;
    synapsePrediction = uint8(3.*ones(length(interfaceSurfaceList),1));
    type1Indices = scores(1:end/2) > 0;
    type2Indices = scores(end/2 + 1:end) > 0;
    synapsePrediction(type1Indices) = 1;
    synapsePrediction(type2Indices) = 2;
catch
    synapsePrediction = zeros(length(interfaceSurfaceList),1,'uint8');
end

data_sa.raw = raw;
data_sa.subsegmentsList = subsegmentsList;
data_sa.segments = segments;
data_sa.neighborIDs = data.neighborIDs;
data_sa.interfaceSurfaceList = interfaceSurfaceList;
data_sa.centroidList = cellfun(@(x)x + [100, 100, 40],data.centroidList,'UniformOutput',false);

end

