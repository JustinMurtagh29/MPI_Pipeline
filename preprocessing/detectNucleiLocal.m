function nuclei = detectNucleiLocal( raw, vessels )
% Detect nuclei on small 3D cubes

% Smooth raw data with gaussian kernel
raw = uint8(smooth3(raw, 'gaussian', 9, 4));

% Calculate gradient magnitude
gradmag = filter3d.gaussiansmoothedgradmagnitude(raw, 5);

% Extract regions that do not have edges
edges = gradmag > 2 | vessels;
edges = bwareaopen(edges, 1e3);
edges = imclose(edges, makeSphere(9));

% Nuclei are regions without edges
nuclei = bwareaopen(imfill(~edges & raw > 110 & raw < 160, 'holes') , 1e6);

end

