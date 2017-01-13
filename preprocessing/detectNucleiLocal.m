function nuclei = detectNucleiLocal( raw, vessels )
% Detect nuclei on small 3D cubes

% Smooth raw data with gaussian kernel (now close to isotrop in physical scale)
raw = smooth3Aniso(raw, [17 17 7], [7.5 7.5 3]);

% Calculate gradient magnitude
gradmag = gradient3Aniso(raw, [7 7 3], [3 3 1.2]);

% Extract regions that do not have edges
edges = gradmag > 2 | vessels;
edges = bwareaopen(edges, 1e3);
edges = imclose(edges, makeSphere(9));

% Nuclei are regions without edges
nuclei = bwareaopen(imfill(~edges & raw > 110 & raw < 160, 'holes') , 1e6);

end
