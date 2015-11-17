function calculateLocalCoMs(p, i, j, k)

    % Load global segmentation (and cut out bboxSmall)
    thisCubeOffset = p.local(i,j,k).bboxSmall(:,1) - [1; 1; 1]; 
    load([p.local(i,j,k).saveFolder 'segGlobal.mat']);
    thisSeg = seg(1-p.tileBorder(1,1):end-p.tileBorder(1,2), ...
        1-p.tileBorder(2,1):end-p.tileBorder(2,2), ...
        1-p.tileBorder(3,1):end-p.tileBorder(3,2));
    clear seg;
    % Remap IDs to make regionprops work better
    bg = thisSeg == 0;
    minObjId = min(thisSeg(~bg));
    thisSeg = thisSeg - (minObjId - 1);
    thisSeg(bg) = 0;
    clear bg;
    % Calculate and save CoMs
    rp = regionprops(thisSeg, 'Centroid');
    localPos = round(cat(1,rp(:).Centroid));
    localPos = localPos(:, [2 1 3]);
    globalPos = bsxfun(@plus, localPos, thisCubeOffset');
    lengthRP = length(rp);
    save([p.local(i,j,k).saveFolder 'comSmall.mat'], 'minObjId', 'lengthRP', 'globalPos');

end

