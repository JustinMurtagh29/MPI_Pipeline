function segmentation = watershedSeg( affZ, r, h, v )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rmpath('/home/mberning/code/KLEE/');
segmentation = cell([length(r) length(h) length(v)]);
tic;
for i=1:length(r)
    if r(i) ~= 0
        [x,y,z] = meshgrid(-r(i):r(i),-r(i):r(i),-r(i):r(i));
        se = (x/r(i)).^2 + (y/r(i)).^2 + (z/r(i)).^2 <= 1;
        % Opening by reconstruction
        affEroded = imerode(affZ, se);
        affRecon = imreconstruct(affEroded, affZ);
        % Closing by reconstruction
        affReconDilated = imdilate(affRecon, se);
        affReconRecon = imreconstruct(imcomplement(affReconDilated), imcomplement(affRecon));
        affReconRecon = imcomplement(affReconRecon);
    else
        affReconRecon = imcomplement(affZ);
    end
    for j=1:length(h)
        % Changes for viewSeg
%         affHmin = imhmin(imcomplement(affReconRecon), h(j), 26);
        % Changes for viewSeg
        affReconRecon = imcomplement(affReconRecon);
        affHmin = imhmin(affReconRecon, h(j), 26);
        bw1 = imregionalmin(affHmin, 26);
        for k=1:length(v)
            if v(k) ~= 0
                bw2 = bwareaopen(bw1, v(k), 26);
            else
                bw2 = bw1;
            end
            affImposedMin = imimposemin(affHmin, bw2, 26);
            segmentation{i,j,k} = watershed(affImposedMin, 26);
            toc
            display([num2str(i) num2str(j) num2str(k)]);
            tic;
        end
    end
end

end

