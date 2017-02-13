function globalSegMat(p, i, j, k)

    % Load offset and segmentation and extract relevant section
    load(p.local(i,j,k).segFile);
    load([p.local(i,j,k).saveFolder 'localToGlobalSegId.mat']);
    
    segNew = zeros(size(seg), 'uint32');
    for l=1:length(localIds)
        segNew(seg == localIds(l)) = globalIds(l);
    end

    clear seg;
    seg = segNew;

    % Write modified seg with globalIDs
    Util.save([p.local(i,j,k).saveFolder 'segGlobal.mat'], seg);

end

