function h = laplacian(siz)
    if siz < 3
        error('Discrete Laplacian with size per dimension < 3 not defined');
    else
        siz = siz + [2 2 2];
        h = zeros(siz);
        h([1 end],:,:) = 1;
        h(:,[1 end],:) = 1;
        h(:,:,[1 end]) = 1;
        h = bwdist(h);
        h([1 end],:,:) = [];
        h(:,[1 end],:) = [];
        h(:,:,[1 end]) = [];
        if rem(siz,2) == 0  
            mid = (siz-2)/2;
            h(mid:mid+1, mid:mid+1, mid:mid+1) = 0;
            h(mid:mid+1, mid:mid+1, mid:mid+1) = repmat(-sum(h(:)),[2 2 2]);
            h = double(h);
        else
            mid = ceil((siz-2)/2);
            h(mid, mid, mid) = 0;
            h(mid, mid, mid) = sum(-h(:));
            h = double(h);
        end
    end
end