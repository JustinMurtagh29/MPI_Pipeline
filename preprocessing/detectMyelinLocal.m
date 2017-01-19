function myelin = detectMyelinLocal( raw )
% Detect myelin on small 3D cubes

% Smooth raw data with gaussian kernel
raw = smooth3Aniso(raw, [21 21 9], [3 3 1.2]);

% Detect things that are dark
darkThings = raw < 94;

% Because filter creates dark areas close to borders,
% Remove detected darkThings there
darkThings([1:3 size(darkThings,1)-2:size(darkThings,1)],:,:) = 0;
darkThings(:,[1:3 size(darkThings,2)-2:size(darkThings,2)],:) = 0;
darkThings(:,:,[1:3 size(darkThings,3)-2:size(darkThings,3)]) = 0;

% Myelin are large dark regions
myelin = bwareaopen(darkThings, 1e5);

% Remove bright objects (smooth regions in apicals mostly) and
% Regions with high average gradients
statsRaw = regionprops(myelin, {'PixelIdxList'});
for i=1:length(statsRaw)
    deciles = quantile(raw(statsRaw(i).PixelIdxList),11);
    if deciles(3) > 85
        myelin(statsRaw(i).PixelIdxList) = 0;
    end
end

end
