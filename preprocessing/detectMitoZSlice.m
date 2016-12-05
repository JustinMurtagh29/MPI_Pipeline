function mito = detectMitoZSlice( raw, borderLoaded )
% Detect myelin on small 3D cubes (typically <=128 z-slices)
% Cuts away borderLoaded

% Smooth raw data with gaussian kernel
raw = uint8(smooth3(raw, 'gaussian', 5, 1.5));

% Detect things that are dark
darkThings = raw < 105;
% Keep only 
mito = bwareaopen(darkThings, 1e4);
mito = imclose(mito, makeSphere(3));
mito = mito(1+borderLoaded(1):end-borderLoaded(1), 1+borderLoaded(2):end-borderLoaded(2), 1+borderLoaded(3):end-borderLoaded(3));

end

function sphere = makeSphere(rad)
    % Based on code by Benedikt Staffler
    % From Util.getPointsInBall
    
    r2 = ceil(rad);
    [xx, yy, zz] = meshgrid(-r2:r2, -r2:r2, -r2:r2);
    sphere = sqrt(xx .^ 2 + yy .^ 2 + zz .^ 2) <= rad;
end
