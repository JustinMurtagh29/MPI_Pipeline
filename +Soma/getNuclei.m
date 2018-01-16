function nuclei = getNuclei(somaIDs, overlapF)
%GETNUCLEI Get the nuclei segment equivalence classes for the L4 dataset.
% INPUT somaIDs: (Optional) [Nx1] int
%           List of the nuclei ids for which the segments are extracted.
%           (Default: all 125 nuclei).
%       overlapF: (Optional) double
%           Fraction of overlap that a segment must have with the nucleus
%           in order to be considered part of the nucleus.
%           (Default: 0.5)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

p = Gaba.getSegParameters('ex145_ROI2017');
meta = load(p.svg.segmentMetaFile, 'voxelCount');

if ~exist('somaIDs', 'var') || isempty(somaIDs)
    somaIDs = 1:125;
end

if ~exist('overlapF', 'var') || isempty(overlapF)
    overlapF = 0.5;
end

% mag stuff to get wk coords
mag1bbox = [128, 128, 128, 5573, 8508, 3413]; 
mag1bbox = Util.convertWebknossosToMatlabBbox(mag1bbox);
mag1bbox = mag1bbox + [-25 25; -25 25; -10 10];

nuclei = cell(length(somaIDs), 1);

Util.log('Loading nuclei masks and segment IDs.');
fprintf('Processing ');
for somaID = somaIDs(:)'
    
    % Note (BS): Load directly from nuclei segmentation?
    m = load(strcat(['/gaba/u/mberning/results/pipeline/20170217_ROI/' ...
        'soma/Nuclei/Nucleus'], int2str(somaID), '.mat'));
    nucleus = m.nucleus;
    if isfield(m, 'bbox') % load bbox from file if exists
        bbox = m.bbox;
    else
        margin = [15, 15, 10]; 
        bboxM4Cropped = [round(rp(somaID).BoundingBox(1:3)) - margin; ...
            round(rp(somaID).BoundingBox(1:3) + ...
            rp(somaID).BoundingBox(4:6)) + margin]';
        bbox = bsxfun(@plus,bsxfun(@plus,4.*bboxM4Cropped,[-1, 2]),  ...
            mag1bbox(:,1));
        % fix of bbox from nuclear_pores/+MouseROI2016/ProcessSingleNucleus2.m
        raw = readKnossosRoi(p.raw.root, p.raw.prefix, bbox);
        if min(min(min(raw))) == 0
            boxsize = size(raw)';
            mask = raw==0;
            maskrp = regionprops(permute(~mask, [2 1 3]));
            maskbox = [ceil(maskrp.BoundingBox(1:3))', ...
                floor(maskrp.BoundingBox(1:3) + maskrp.BoundingBox(4:6))'];
            bbox(:,1) = bsxfun(@plus, ...
                bsxfun(@plus, bbox(:,1), maskbox(:,1)), -1);
            bbox(:,2) = bsxfun(@plus, bbox(:,2), -(boxsize - maskbox(:,2)));
        end
        clear raw mask maskrp maskbox boxsize
    end
    
    warning('off', 'all');
    somaSeg = Seg.IO.loadSeg(p, bbox);
    warning('on', 'all');

    % this caused some bugs so assert it
    assert(all(size(somaSeg) == size(nucleus)));
    
    % get nucleus seg ids with at least 50% overlap with nucleus
    nucleusSegIds = tabulate([0; somaSeg(nucleus)]); % add 0 just to make sure
    nucleusSegIds = nucleusSegIds(2:end, 1:2);
    toKeep = nucleusSegIds(:,2) > overlapF*meta.voxelCount(nucleusSegIds(:,1));
    nucleusSegIds = nucleusSegIds(toKeep, 1);
    nuclei{somaID} = nucleusSegIds;
    clear nucleus somaSeg
    fprintf('.');
end
fprintf(' done\n');
Util.log('Finished loading nuclei.');

end