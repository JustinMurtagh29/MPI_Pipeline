function globalSegId(p, coords, i, j, k)

    % Load offset and segmentation and extract relevant section
    load(p.local(i,j,k).segFile);
    load([p.saveFolder 'numEl.mat']);
    seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));
    
    % Linearize segmentation for relabeling
    sizeCube = size(seg);
    seg = seg(:);

    % Find mapping from (non-continous due to removal of overlap) local IDs to continous global IDs
    [localIds, ~, n] = unique(seg);
    
    % Take into account special role of 0 as background ID which remains unchanged
    nrUniqueIds = length(localIds) - 1;
    if any(localIds == 0)
        globalIds = [uint32(0) uint32(1:nrUniqueIds) + numElTotal(i,j,k)];
    else
        % In case of only one segment, otherwise there should be borders, throw error if this is not the case
        if nrUniqueIds ~= 0
           error('Something wrong');
        else
            globalIds = uint32(1) + numElTotal(i,j,k);
        end 
    end

    % This is the array which saves the mapping from local to global IDs
    localToGlobalMapping = [localIds globalIds'];

    % Add sum of all former 
    seg = globalIds(n);
    seg = reshape(seg, sizeCube);

    % Write modified seg with globalIDs (not complete wrt IDs, e.g. localIDs in outer bounidng box of inner cube will not be present, any better idea?)
    writeKnossosRoi(p.seg.root, p.seg.prefix , [p.local(i,j,k).bboxBig(:,1) + coords{i,j,k}(:,1) - [1; 1; 1]]', seg, 'uint32'); 
    % Save mapping from local to global IDs
    save([p.local(i,j,k).saveFolder 'localToGlobalSegId.mat'], 'localToGlobalMapping');

end

