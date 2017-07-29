function subsegmentsList = calculateSubsegments( contactSurfaceList,neighborIDs,segments,rinclude )
%CALCULATESUBSEGMENTS Calculate subsegments in parallel.

fprintf('%s INFO restricting segments to parts close to contact\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));

segmentIndices = SynEM.Paper.TrainingSet.bucketsort(segments(:) + 1);
segmentIndices = segmentIndices(2:end);


subsegmentsList = cell(1,length(rinclude));
for k = 1:length(subsegmentsList)
    subsegmentsList{k} = cell(length(contactSurfaceList),2);
end

for iContact = 1:length(contactSurfaceList)
    %get contact wall indices in nm
    [x,y,z] = ind2sub(size(segments),contactSurfaceList{iContact});
    wallIndicesNM = bsxfun(@times,[x,y,z] - 1,[11.24,11.24,28]);

    % every location out of this from cannot be within rinclude to any contact voxel
    roughframe = [min(wallIndicesNM,[],1)-max(rinclude); max(wallIndicesNM,[],1)+max(rinclude)];

    %slice neighborIDs
    currentNeighbors = neighborIDs(iContact,:);

    for iCell = 1:2
        %calculate segment indices in nm
        currentSegmentIndices = segmentIndices{currentNeighbors(iCell)};
        [x,y,z] = ind2sub(size(segments),currentSegmentIndices);
        segmentIndicesNM= bsxfun(@times,[x,y,z] - 1,[11.24,11.24,28]);

        %discard voxels outside the roughframe
        sortedout = any(bsxfun(@lt, segmentIndicesNM, roughframe(1,:)), 2)|any(bsxfun(@gt, segmentIndicesNM, roughframe(2,:)), 2); 
        currentSegmentIndices = currentSegmentIndices(~sortedout);
        segmentIndicesNM = segmentIndicesNM(~sortedout,:);

        %calculate subsegments using rangesearch
        [idx,D] = rangesearch(segmentIndicesNM,wallIndicesNM,max(rinclude));
        D = cell2mat(D');
        idx = cell2mat(idx');
        for k = 1:length(rinclude)
            indices = D < rinclude(k);
            %quicker unique
            indices = sort(idx(indices));
            indices = indices([true;diff(indices(:))>0]);
            
            subsegmentsList{k}{iContact,iCell} = currentSegmentIndices(indices,:);
        end
    end
end

end

