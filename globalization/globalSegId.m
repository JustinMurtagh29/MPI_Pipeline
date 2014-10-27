function globalSegId(p, coords, i, j, k)

    % Load offset and segmentation and extract relevant section
    load(p.local(i,j,k).segFile);
    load([p.saveFolder 'numEl.mat']);
    seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));

    % Find mapping from (non-continous due to removal of overlap) local IDs to continous global IDs
    seg = seg(:);
    [localIds, ~, n] = unique(seg);
    globalIds = 1:length(localIds) + numElTotal(i,j,k);

    % This is the array which saves the mapping from local to global IDs
    localToGlobalMapping = [localIds; globalIds];

    % Add sum of all former 
    seg = globalIds(;
    seg = reshape(seg, [coords{i,j,k}(:,2)-coords{i,j,k}(:,1) + 1]');

    % Write modified seg with globalIDs (not complete wrt IDs, e.g. localIDs in outer bounidng box of inner cube will not be present, any better idea?)
    writeKnossosRoi(p.seg.root, p.seg.prefix , [p.local(i,j,k).bboxBig(:,1) + coords{i,j,k}(:,1) - [1; 1; 1]]', seg, 'uint32'); 
    % Save mapping from local to global IDs
    save([p.local(i,j,k).saveFolder 'localToGlobalSegId.mat'], 'localToGlobalMapping');

end

