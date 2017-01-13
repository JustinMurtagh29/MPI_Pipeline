function myelin = detectMyelinLocal( raw )
% Detect myelin on small 3D cubes

% Smooth raw data with gaussian kernel
raw = smooth3(raw, [7 7 3], [3 3 1.2]);

% Detect things that are dark
darkThings = raw < 93;

% Myelin are large dark regions
myelin = bwareaopen(darkThings, 1e5);

end


