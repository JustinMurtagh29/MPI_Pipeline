function [contactSurfaceList, subsegmentsList,neighborIDs] = contactCalculationPipeline(segments, walls, rinclude,bbox,useInGUI)
% CONTACTCALCULATIONPIPELINE Generates lists of indices of voxels from the contact surfaces and close pre and postsynaptic parts give a segmentation
% 
% INPUT segments:Image in which all segments have an unique integer ID
%       walls: Binary image of the same size as segments indicating walls
%              with true
%       rinclude: maximum distance in nm a voxel can have to be included in
%                 the pre and post synaptic segment parts (can be a vector
%                 of multiple distances)
%       bbox: bounding box of raw
% 
% OUTPUT contactSurfaceList: Cell array of size n by 1, where n is the number of 
%                    contacts. Each element is m_i by 1 array in which each
%                    element represents a linear index of a voxel of
%                    contact i
%        subsegmentsList: Cell array with length(rinclude) entries each
%                    containing a n times 2 cell array conaining the linear
%                    indices of voxels of a segment restricted to a radius
%                    of rinclude within the contact surface
%        neighborIDs: double array containing the IDs of the neighboring
%                     segments of the correspoding contact
% 
% created by: Jan Gleixner (jan.gleixner@gmail.com)
% modified by: Benedikt Staffler

fprintf('%s INFO starting getIdx\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));
 
[ wallIndices, uneighboursID] = SynEM.Paper.TrainingSet.getContactIndices(segments, walls); 
if isempty(wallIndices)
    error('No wall voxels could be found. Aborting getContacts.');
end
 
% generate list of wall indxi of contacts
% first: sorting all wall indices from contacts according to cell 1 and split
left = uneighboursID(:,2);
right = uneighboursID(:,1);
[sleft, left2sleft] = sort(left); % should be bucket sort
usleftidxl = [true; diff(sleft)>0]; % first of occurrences of unique ids
usleft = sleft(usleftidxl); % get unique ids
usleftidxi = find(usleftidxl); % idxi of first of occurrences of unique ids
leftSortedWallIndices = wallIndices(left2sleft); % sort wallvoxelidxi in the same way as left
leftSortedWallIndices = mat2cell(leftSortedWallIndices, ...
    diff([usleftidxi; numel(usleftidxl)+1]), 1); % and separate
leftSortedRight = right(left2sleft); % sort rightids in the same way as left
leftSortedRight = mat2cell(leftSortedRight, ...
    diff([usleftidxi; numel(usleftidxl)+1]), 1); % and separate

 
% second: sorting all wall idxi from all contacts of each cell 1
% according to cell 2 and split
% [sright, right2sright] = cellfunn(2, @(lsright)sort(lsright), lsright);
[sright, right2sright] = cellfun(@(lsright)sort(lsright), ...
    leftSortedRight,'UniformOutput',false);
usrightidxl = cellfun(@(sright)[true; diff(sright)>0], sright, ...
    'UniformOutput', false); % first of occurrences of unique ids
usright = cellfun(@(sright, usrightidxl)sright(usrightidxl), ...
    sright, usrightidxl, 'UniformOutput', false); % get unique ids
usrightidxi = cellfun(@(usrightidxl)find(usrightidxl), usrightidxl, ...
    'UniformOutput', false); % idxi of first of occurrences of unique ids
contactSurfaceList = cellfun(@(lswallidxi, right2sright) ...
    lswallidxi(right2sright), leftSortedWallIndices, right2sright, ...
    'UniformOutput', false); % sort wallvoxelidxi in the same way as right
contactSurfaceList = cellfun(@(rlswallidxi, usrightidxi, usrightidxl) ...
    mat2cell(rlswallidxi, diff([usrightidxi; numel(usrightidxl)+1]), 1),...
    contactSurfaceList, usrightidxi, usrightidxl, 'UniformOutput', false); % and separate
contactSurfaceList = cat(1,contactSurfaceList{:});
leftright = arrayfun(@(i)[repmat(usleft(i), size(usright{i})), usright{i}], ...
    1:numel(usright), 'UniformOutput', false); 
leftright = cat(1, leftright{:});

fprintf('%s INFO found %s contacts\n', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'), ...
    mat2str(length(contactSurfaceList)));

[contactSurfaceList,neighborIDs] = SynEM.Paper.TrainingSet.splitDiscontiguousContacts(contactSurfaceList,leftright,size(segments),bbox,useInGUI);

fprintf('%s INFO found %s contact surfaces\n', ...
    datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'), ...
    mat2str(length(contactSurfaceList)));

subsegmentsList = SynEM.Paper.TrainingSet.calculateSubsegments( ...
    contactSurfaceList,neighborIDs,segments,rinclude);
end