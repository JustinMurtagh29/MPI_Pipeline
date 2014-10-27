function globalSegId(p, i,j,k)

    coords = localToGlobalBBoxModification(p);
    load([p.seg.root 'numEl.mat'])

    %Load segmentation and extract relevant section
    load(p.local(i,j,k).segFile)
    seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));

    % 
    [uni, ia, ic] = unique(seg(:));
    nrGlobalIDs = length(uni)-1;
    % give new segmentIds from 1:N. Minus 1 bc. of zero segment
    segNew = ic - 1;

    clear UV1 UV2
    %find unique values without rearrangment. No 'sorted' option in Matlab 2011
    [X1, sortVec1] = sort(seg(:));
    [X2, sortVec2] = sort(segNew);
    UV1(sortVec1) = ([1; diff(X1)] ~= 0); 
    UV2(sortVec2) = ([1; diff(X2)] ~= 0);
    Y1 = seg(UV1);
    Y2 = segNew(UV2);

    %this is the globalSegId Array which gives every 'old' segment it's new globalID
    checkMat = [Y1', Y2];

    % Add sum of all former 
    seg = uint32(segNew) + numElTotal(i,j,k);
    checkMat = uint32(checkMat) + numElTotal(i,j,k);
    seg = reshape(seg, [coords{i,j,k}(:,2)-coords{i,j,k}(:,1) + 1]');

    % Write modified seg with globalIDs (not complete wrt IDs, e.g. localIDs in outer bounidng box of inner cube will not be present, any better idea?)
    writeKnossosRoi(p.seg.root, p.seg.prefix , [p.local(i,j,k).bboxBig(:,1) + coords{i,j,k}(:,1) - [1; 1; 1]]', seg, 'uint32'); 
    save([p.local(i,j,k).saveFolder 'globalSegId.mat'], 'checkMat');

end

