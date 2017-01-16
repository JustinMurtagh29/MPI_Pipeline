function nuclei = detectNucleiLocal( raw, vessels )
% Detect nuclei on small 3D cubes

% Smooth raw data with gaussian kernel (now close to isotrop in physical scale)
raw = smooth3Aniso(raw, [21 21 9], [3 3 1.2]);

% Calculate gradient magnitude
gradmag = gradient3Aniso(raw, [11 11 5], [2 2 0.8]);

% Extract regions that have edges
edges = gradmag > 5.5 | vessels;
edges = bwareaopen(edges, 1e3);
edges = imclose(edges, makeSphere([10 10 4], 10));

% Nuclei are regions without edges
nuclei = bwareaopen(imfill(~edges & raw > 120 & raw < 150, 'holes') , 1e6);

end
