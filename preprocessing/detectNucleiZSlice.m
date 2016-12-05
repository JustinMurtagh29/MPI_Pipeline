function nuclei = detectNucleiZSlice( raw, borderLoaded )
% Detect nuclei on small 3D cubes (typically <=128 z-slices)
% Cuts away borderLoaded

% Smooth raw data with gaussian kernel
raw = uint8(smooth3(raw, 'gaussian', 11, 4));

% Calculate gradient magnitude
gradmag = filter3d.gaussiansmoothedgradmagnitude(raw, 5);
% Remove 1st part border loaded in addition (keep rest for larger objects
% during morphological operations)
toRemoveFirst = [5 5 5];
gradmag = gradmag(1+toRemoveFirst(1):end-toRemoveFirst(1), 1+toRemoveFirst(2):end-toRemoveFirst(2), 1+toRemoveFirst(3):end-toRemoveFirst(3));

% Search for 2 regions, 1 with nuclei one without and look at histogram
%a = gradmag(1:100, 1:100, 1:28); b = gradmag(end-99:end, end-99:end, 1:28);
%figure; subplot(2,1,1); histogram(a(:), 0:0.5:10); subplot(2,1,2); histogram(b(:), 0:0.5:10);

% Extract regions that do not have (many continouus) edges
edges = gradmag > 2;
edges = bwareaopen(edges, 1000);
edges = imclose(edges, makeSphere(10));

%
nuclei = bwareaopen(~edges, 1e6);
toRemoveSecond = borderLoaded - toRemoveFirst;
nuclei = nuclei(1+toRemoveSecond(1):end-toRemoveSecond(1), 1+toRemoveSecond(2):end-toRemoveSecond(2), 1+toRemoveSecond(3):end-toRemoveSecond(3));

end

function sphere = makeSphere(rad)
    % Based on code by Benedikt Staffler
    % From Util.getPointsInBall
    
    r2 = ceil(rad);
    [xx, yy, zz] = meshgrid(-r2:r2, -r2:r2, -r2:r2);
    sphere = sqrt(xx .^ 2 + yy .^ 2 + zz .^ 2) <= rad;
end
