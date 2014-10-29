function segmentation = watershedSeg_v3( aff, thres, v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rmpath('/home/mberning/code/KLEE/');
segmentation = cell([length(thres) length(v)]);
tic;
aff = imcomplement(aff);
for i=1:length(thres)
    fgm1 = affSmooth < thres(i);
    fgm2 = imfill(fgm1,'holes');
    for k=1:length(v)
        fgm3 = bwareaopen(fgm2, v(k), 26);
        affMarker = aff;
        affMarker = imimposemin(affMarker, fgm3);
        segmentation{i,k} = watershed(affMarker, 26);
        display([num2str(i) num2str(k) ':']);
        toc
        tic;
    end
end
toc

end

