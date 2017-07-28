function [rp, fpCen, fnCen, confSeg, fnIds] = segmentationOverlapRestricted( ...
    synSeg, pred, synIds, border, thresholds, sizeTLower, sizeTUpper, ...
    overlapT, mergeTerm)
%SEGMENTATIONOVERLAPRESTRICTED Calculate prediction performance by counting
%tp/fn detections for the specified synapses only and fps only for
%predictions in a smaller bounding box.
% INPUT synSeg: n-d int
%           The ground truth segmentation with ids for each synapse.
%       pred: n-d float
%           The voxel-wise prediction from which proposal segmentations are
%           calculated.
%       synIds: [Nx1] int
%           The ids in synSeg that indicate objects that need to be
%           found. Ids that are not in synIds but in synSeg will not
%           cause detections errors. This was done on purpose to allow for
%           a kind of masking.
%           E.g. to only check for excitatory synapses in synSeg specify
%           the ids of the excitatory synapses in synIds but keep the
%           inhibitory synapse ids in synSeg. Similarly, to minimize
%           boundary effects of a volume this can be used to select only
%           synapses in the center of a bounding box while synapses at the
%           border won't cause detection errors.
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
%        fpCen: [Nx3] int
%           FP detection centroids w.r.t. synSeg.
%           (only applicable for a single threshold).
%        fnCen: [Nx3] int
%           FN detection centroids w.r.t. synSeg.
%           (only applicable for a single threshold).
%        confSeg: 3d uint8
%           Label matrix containing confusion matrix labels for each
%           predicted component and FN component. Labels are
%           '1': TP
%           '2': FP
%           '3': FN
%           (only applicable for a single threshold)
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
% if nargout > 1 && length(thresholds) > 1
%     error('Output arguments 2-4 can only be calculated for a single threshold.');
% end
if nargout > 1
    confSeg = zeros(size(pred), 'uint8');
end
    

rp = zeros(length(thresholds), 2);
fnIds = cell(length(thresholds), 1);
tic;
for i = 1:length(thresholds)
    % get predicted segmentation
    cPred = SynEM.Eval.postProcessing(pred, thresholds(i), ...
        sizeTLower, sizeTUpper);
    cPredL = bwlabeln(cPred);
    
    stats = regionprops(cPredL, 'PixelIdxList');
    
    % get synapse indices for each pred component
    predGTOverlap = arrayfun(@(x)getOverlapp(synSeg, ...
        x.PixelIdxList, overlapT), ...
        stats, 'UniformOutput', false);
    synOverlapIdx = ~cellfun(@isempty, predGTOverlap);
    
    % delete synaptic components (such that only FPs remain)
    cPredL(ismember(cPredL, find(synOverlapIdx))) = 0;
        
    if mergeTerm && any(cellfun(@length, predGTOverlap) > 1) % terminate
        rp(i,:) = NaN;
        continue;
    end
    
    % tps & fps
    fns = setdiff(synIds, cell2mat(predGTOverlap(synOverlapIdx)));
    tps = length(setdiff(synIds, fns));
    if nargout > 1
        confSeg(ismember(synSeg, fns)) = 3; % fns
        confSeg(ismember(synSeg, setdiff(synIds, fns))) = 1; % fps
    end
    fnIds{i} = fns;
    fns = length(fns);
    
    % fps:
    %restricting pred to the smaller cube
    if size(border, 1) > 1
        toDel = false(size(cPredL) - sum(border, 1));
        toDel = padarray(toDel, border(1,:), true, 'pre');
        toDel = padarray(toDel, border(2,:), true, 'post');
        cPredL(toDel) = 0;
    else
        cPredL(padarray(false(size(cPredL) - 2.*border), border, true)) = 0;
    end
    
    % delete small components
    cPredL = bwlabeln(cPredL > 0); %renumber starting from one
    if nargout > 1
        stats = regionprops(cPredL, 'Area', 'PixelIdxList', 'Centroid');
    else
        stats = regionprops(cPredL, 'Area');
    end
    
    % all remaining components are fps (snypases were deleted above)
    fps = ([stats.Area] > sizeTLower) & ([stats.Area] < sizeTUpper);
    if nargout > 1
        confSeg(cell2mat({stats(fps).PixelIdxList}')) = 2; % fps
        fpCen = round(reshape([stats(fps).Centroid], 3, [])');
        fpCen = fpCen(:,[2 1 3]);
    end
    fps = sum(fps);
             
    % calculate rp
    rp(i,1) = tps/(tps + fns); %recal
    rp(i,2) = tps/(tps + fps); %precision
    Util.progressBar(i, length(thresholds));
end

if nargout > 1
    stats = regionprops(confSeg == 3, 'Centroid');
    fnCen = round(reshape([stats.Centroid], 3, [])');
    fnCen = fnCen(:,[2 1 3]);
end
end

function synIds = getOverlapp(synSeg, pixelIdxList, overlapT)
%Get the ground truth IDs at the pixelIdsList locations.
ids = synSeg(pixelIdxList);
%get only ids with at least overlapT occurrences
synIds = find(accumarray(ids(ids > 0), 1) > overlapT);

end

