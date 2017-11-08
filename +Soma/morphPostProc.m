function somaAgglos = morphPostProc( p, somaAgglos )
%MORPHPOSTPROC Morphological post-processing of the somas to close holes
% and get rid of extensive dendrite exits.
% INPUT p: struct
%           Segmentation parameter struct.
%       somaAgglos: [Nx1] cell
%           Cell array of soma agglos. Each cell contains a list of integer
%           id of segments belonging to a single soma.
% OUTPUT somaAgglos: [Nx1] cell
%           The updated soma agglos.
%
% NOTE This function requires ~60G of RAM on the ex145 dataset.
%
% Based on code by: Kevin Boergens <kevin.boergens@brain.mpg.de>
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% load mag8 segmentation of whole dataset
dat = Datasets.WkDataset(p.seg);
dat8 = Datasets.WkDataset(p.seg);
dat8 = dat8.switchMag(8);
bboxSegMag8 = Datasets.WkDataset.transformCoordsToMag(p.bbox', 1, 8)';
fS = [4, 4, 2]; % choose bbox size to be divisible by filter size
bboxSegMag8(:,2) = bboxSegMag8(:,2) + (fS(:) - rem(diff(bboxSegMag8, [], 2) + 1, fS(:)));
segMag8 = dat8.readRoi(bboxSegMag8);

numSomas = length(somaAgglos);
for idx = 1 : numSomas
    
    % morphological operations to close holes and get rid of exits
    somaMaskMag8 = ismember(segMag8, somaAgglos{idx});
    somaMaskMag8 = nlfilter3(somaMaskMag8, @max, fS);
    somaMaskMag8 = imopen(somaMaskMag8, ones([7,7,7]));
    somaMaskMag8 = imdilate(somaMaskMag8, ones([3,3,3]));
    
    % upsample mask 8 to mag 8 again
    somaMaskMag8 = Soma.resize(somaMaskMag8, size(segMag8));
    
    % get cube in mag1
    stats = regionprops(somaMaskMag8, 'BoundingBox', 'Area');
    [~, toKeep] = max([stats.Area]); % keep only largest component
    stats = stats(toKeep);
    bboxMag8 = reshape(stats.BoundingBox, [], 2);
    bboxMag8(:,1) = ceil(bboxMag8(:,1));
    bboxMag8(:,2) = bboxMag8(:,1) + bboxMag8(:,2) - 1;
    bboxMag8 = bboxMag8([2 1 3],:);
    
    % load segmentation in bbox
    bboxMag1 = Datasets.WkDataset.transformCoordsToMag(bboxMag8', 8, 1);
    bboxMag1 = bboxMag1';
    bboxMag1 = bboxMag1 + 128;
    seg = dat.readRoi(bboxMag1);
    somaMaskMag8 = Util.getBbox(somaMaskMag8, bboxMag8);
    mask = Soma.resize(somaMaskMag8, diff(bboxMag1, [], 2)' + 1);
    mask = imerode(mask, ones(3, 3, 3)); % just to make sure it is not too large
    somaAgglos{idx} = setdiff(seg(mask), 0);
end

end
