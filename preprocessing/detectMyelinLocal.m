function myelin = detectMyelinLocal( raw )
% Detect myelin on small 3D cubes

% Smooth raw data with gaussian kernel
raw = uint8(smooth3(raw, 'gaussian', 5, 3));

% Detect things that are dark
darkThings = raw < 95;

% Myelin are large dark regions
myelin = bwareaopen(darkThings, 1e5);

end


