function [somaAgglos, toDiscardSomaSegs, toDiscardSyn] = ...
    morphPostProc( p, somaAgglos, segmentCom, synComs )
%MORPHPOSTPROC Morphological post-processing of the somas to get rid of
% extensive dendrite exits.
% INPUT p: struct
%           Segmentation parameter struct.
%       somaAgglos: [Nx1] cell
%           Cell array of soma agglos. Each cell contains a list of integer
%           id of segments belonging to a single soma.
%       segmentCom: [Nx3] double
%           Global segment com list. This is used to keep only soma
%           segments whose com is within the constructed soma mask.
%       synComs: (Optional) [Nx1] cell of [Mx3] double
%           Cell array with synapse centroids for each synapse. Synapse
%           centroids are filtered as segment coms if they are provided.
% OUTPUT somaAgglos: [Nx1] cell
%           The updated soma agglos.
%        toDiscardSomaSegs: [Nx1] cell of [Mx1] logical
%           Logical list of soma segments for each soma that are not in the
%           morphological mask and could potentially be discarded.
%        toDiscardSyn: [Nx1] cell of [Mx1] logical
%           Logical list of synapses for each soma that are not in the
%           morphological mask and could potentially be discarded.
% Based on code by: Kevin Boergens <kevin.boergens@brain.mpg.de>
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% load mag8 segmentation of whole dataset
dat8 = Datasets.WkDataset(p.seg);
dat8.switchMag(8);
bboxSegMag8 = Datasets.WkDataset.transformCoordsToMag(p.bbox', 1, 8)';
bboxSegMag8(:,1) = 1;
fS = [4, 4, 2]; % choose bbox size to be divisible by filter size
bboxSegMag8(:,2) = bboxSegMag8(:,2) + ...
    (fS(:) - rem(diff(bboxSegMag8, [], 2) + 1, fS(:)));
segMag8 = dat8.readRoi(bboxSegMag8);
coordConvFact = 1./(8.*fS);

numSomas = length(somaAgglos);
toDiscardSomaSegs = cell(numSomas, 1);
toDiscardSyn = cell(numSomas, 1);
for idx = 1 : numSomas
    idx
    % morphological operations to close holes and get rid of exits
    somaMaskMag8 = ismember(segMag8, somaAgglos{idx});
    somaMaskMag8 = nlfilter3(somaMaskMag8, @max, fS);
    somaMaskMag8 = imopen(somaMaskMag8, ones([7,7,7]));
    somaMaskMag8 = imdilate(somaMaskMag8, ones([3,3,3]));
    
    % keep only soma segments that are within the mask
    thisSegComsMag8 = Util.sub2ind(size(somaMaskMag8), ...
        round(bsxfun(@times, segmentCom(somaAgglos{idx},:), coordConvFact)));
    comsInMask = somaMaskMag8(thisSegComsMag8);
    toDiscardSomaSegs{idx} = ~comsInMask;
    somaAgglos{idx}(~comsInMask) = [];
    
    % same for synapse coms
    if exist('synComs', 'var') && ~isempty(synComs) && ~isempty(synComs{idx})
        thisSynComsMag8 = Util.sub2ind(size(somaMaskMag8), ...
            round(bsxfun(@times, synComs{idx}, coordConvFact)));
        comsInMask = somaMaskMag8(thisSynComsMag8);
        toDiscardSyn{idx} = ~comsInMask;
    end
end

end
