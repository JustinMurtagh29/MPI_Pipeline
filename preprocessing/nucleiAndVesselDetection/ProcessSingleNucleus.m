function ProcessSingleNucleus(bbox, mask, vessel, root, prefix, id)
% Processing of a single nucleus in magnification1 in order to detect a
% mask describing the extension of the nucleus.
seg_root = strrep(root,'corrected/color','nuclei/segmentation');
seg_prefix = strrep(prefix,'corrected','nuclei');

gradient_threshold = 0.025;


    raw = readKnossosRoi(root, prefix, bbox);
    %% Elimination of dark regions originating from bleaching
    % Due to the gradient detection which is used to detect the nuclei the
    % borders of these dark regions have to be eliminated.
    randoms=randi([170 205],size(raw,1),size(raw,2),size(raw,3));
    darkregion = raw < 100;
    darkregion = imclose(darkregion,strel('disk',4,0));
    for k=1:size(darkregion,3)
        darkregion(:,:,k) = bwareaopen(darkregion(:,:,k), 100);
    end
    darkregion = imdilate(darkregion,strel('disk',5,0));
    darkregion = uint8(smooth3(darkregion, 'gaussian', 7, 3));
    
    darkborder = darkregion;
    darkregion = logical(darkregion);
        
    raw(darkregion) = randoms(darkregion);
    clear darkregion
    
    [Gx, Gy, Gz] = NuclearPores.imgradientxyz(darkborder);
    darkborder = Gx.^2 + Gy.^2 + Gz.^2;
    clear Gx Gy Gz
    maximum = max(max(max(darkborder)));
    darkborder = darkborder./maximum;
    darkborder = darkborder > gradient_threshold;
    
    darkborder = bwareaopen(darkborder,4000,8);
    darkborder = imdilate(darkborder, strel('disk',3,0));
    
    %% Gradient detection
    raw = uint8(smooth3(raw, 'gaussian', 5, 2));
    [Gx, Gy, Gz] = NuclearPores.imgradientxyz(raw);
    gradsum = Gx.^2 + Gy.^2 + Gz.^2;
    clear Gx Gy Gz
    maximum = max(max(max(gradsum)));
    gradsum = gradsum./maximum;
    %%
    nucleus = gradsum > gradient_threshold;
    clear gradsum
    % Elimination of the edges originating from the dark regions
    nucleus(darkborder) = 0;
    for k=1:size(nucleus,3)
        nucleus(:,:,k) = bwareaopen(nucleus(:,:,k), 250);
    end
    nucleus = bwareaopen(nucleus, 1000);
    clear darkborder
    %% Elimination of mask and vessel segments
    mask = NuclearPores.resize(mask, size(mask).*[4 4 4]);
    nucleus(mask) = 1;
    vessel = NuclearPores.resize(vessel, size(vessel).*[4 4 4]);
    nucleus(vessel) = 1;
    clear vessel
    %%
    nucleus = imclose(nucleus, NuclearPores.makeSphere(10));
    nucleus = bwareaopen(~nucleus, round(0.01*(sum(sum(sum(nucleus))))));
    %%
    nucleus = imclose(nucleus,strel('disk',13,0));
    
    for k=1:size(nucleus,3)
        nucleus(:,:,k) = ~bwareaopen(~nucleus(:,:,k),3000,4);
    end
    
    nucleusMasked = or(nucleus,mask);
    nucleus = imfill(nucleusMasked, 'holes') & ~mask;
    clear mask
    nucleus = bwareaopen(nucleus, round(0.3 * sum(sum(sum(nucleus)))));
    
    directory = strcat('/gaba/u/scchr/MouseP28/NucleiVesselDetection/singleMag1/Nucleus',int2str(id), '.mat');
    save(directory, 'nucleus', '-v7.3')