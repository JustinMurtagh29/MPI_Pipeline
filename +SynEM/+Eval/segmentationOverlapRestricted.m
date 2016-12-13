function rp = segmentationOverlapRestricted( synSeg, pred, ...
    synIds, border, thresholds, sizeTLower, sizeTUpper, overlapT, mergeTerm)
%SEGMENTATIONOVERLAPRESTRICTED Calculate prediction performance by counting
%tp/fn detections for the specified synapses only and fps only for
%predictions in a smaller bounding box.
% INPUT synSeg: n-d int
%           The ground truth segmentation with ids for each synapse.
%       pred: n-d float
%           The voxel-wise prediction from which proposal segmentations are
%           calculated.
%       border: [1x3] int
%           The size of between inner and outer synSeg cube. FPs will only
%           be counted in the inner synSeg cube.
%       thresholds: [Nx1] double
%           List of thresholds for pred to produce the output
%           segmentations.
%       sizeTLower: (Optional) double
%           Lower threshold on the connected component size of the proposal
%           segmentation.
%           (Default: 0)
%       sizeTUpper: (Optional) double
%           Upper threshold on the connected component size of the proposal
%           segmentation.
%           (Default: 0)
%       overlapT: (Optional) double
%           Minimal number of voxels in the overlap between proposal
%           components and ground truth components to accept a tp
%           detection.
%           (Default: 0)
%       mergeTerm: (Optional) logical
%           Flag indicating if that only proposal segmentations are
%           accepted that do not overlap with multiple ground truth
%           components (merge terminate).
%           (Default: true)
% OUTPUT rp: [Nx2] double
%           Recall-precision pairs for each input threshold.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('sizeTLower','var') || isempty(sizeTLower)
    sizeTLower = 0;
end
if ~exist('sizeTUpper','var') || isempty(sizeTUpper)
    sizeTUpper = Inf;
end
if ~exist('overlapT','var') || isempty(overlapT)
    overlapT = 0;
end
if ~exist('mergeTerm','var') || isempty(mergeTerm)
    mergeTerm = true;
end

rp = zeros(length(thresholds), 2);
for i = 1:length(thresholds)
    %get voxel prediction
    cPred = SynEM.Eval.postProcessing(pred, thresholds(i), ...
        sizeTLower, sizeTUpper);
    cPredL = bwlabeln(cPred);
    
    stats = regionprops(cPredL, 'PixelIdxList');
    
    %get synapse indices for each pred component
    predGTOverlap = arrayfun(@(x)getOverlapp(synSeg, ...
        x.PixelIdxList, overlapT), ...
        stats, 'UniformOutput', false);
    synOverlapIdx = ~cellfun(@isempty, predGTOverlap);
    
    %delete synaptic components
    cPredL(ismember(cPredL, find(synOverlapIdx))) = 0;
        
    if mergeTerm && any(cellfun(@length, predGTOverlap) > 1) % terminate
        rp(i,:) = NaN;
        continue;
    end
    
    %tps & fps
    fns = setdiff(synIds, cell2mat(predGTOverlap(synOverlapIdx)));
    tps = length(setdiff(synIds, fns));
    fns = length(fns);
    
    %fps:
    %restricting pred to the smaller cube
    cPredL(padarray(false(size(cPredL) - 2.*border), border, true)) = 0;
    
    %delete small components
    cPredL = bwlabeln(cPredL > 0); %renumber starting from one
    stats = regionprops(cPredL, 'Area');
    
    %all remaining components are fps
    fps = ([stats.Area] > sizeTLower) & ([stats.Area] < sizeTUpper);
    fps = sum(fps);
             
    %calculate rp
    rp(i,1) = tps/(tps + fns); %recal
    rp(i,2) = tps/(tps + fps); %precision
    
end
end

function synIds = getOverlapp(synSeg, pixelIdxList, overlapT)
%Get the ground truth IDs at the pixelIdsList locations.
ids = synSeg(pixelIdxList);
%get only ids with at least overlapT occurrences
synIds = find(accumarray(ids(ids > 0), 1) > overlapT);

end

