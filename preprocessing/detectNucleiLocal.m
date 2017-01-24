function nuclei = detectNucleiLocal( raw, vessels )
% Detect nuclei on small 3D cubes

% Smooth raw data with gaussian kernel (now close to isotrop in physical scale)
raw = smooth3Aniso(raw, [21 21 9], [3 3 1.2]);

% Calculate gradient magnitude
gradmag = gradient3Aniso(raw, [21 21 9], [3 3 1.2]);

% Extract regions that have edges
edges = gradmag > 10 | vessels;
edges = bwareaopen(edges, 1e3);
edges = imclose(edges, makeSphere([10 10 4], 10));

% Nuclei are regions without edges
nuclei = bwareaopen(imfill(~edges, 'holes') , 5e5);

% Remove bright objects (smooth regions in apicals mostly) and
% Regions with high average gradients
statsRaw = regionprops(nuclei, raw, {'Area', 'BoundingBox', 'Centroid', 'MinIntensity', 'MeanIntensity', 'MaxIntensity' 'PixelIdxList'});
statsGradmag = regionprops(nuclei, gradmag, {'Area', 'BoundingBox', 'Centroid', 'MinIntensity', 'MeanIntensity', 'MaxIntensity' 'PixelIdxList'});
for i=1:length(statsRaw)
    if statsRaw(i).MeanIntensity > 150 || statsGradmag(i).MeanIntensity > 6
        nuclei(statsRaw(i).PixelIdxList) = 0;
    end
end

% Remove small edges in detected nuclei
nuclei = imfill(imclose(nuclei, makeSphere([10 10 4], 10)), 'holes');

end

