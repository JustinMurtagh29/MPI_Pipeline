function myelin = detectMyelinLocal( raw )
% Detect myelin on small 3D cubes

% Smooth raw data with gaussian kernel
raw = smooth3Aniso(raw, [11 11 5], [3 3 1.2]);

% Detect things that are dark
darkThings = raw < 92;

% Myelin are large dark regions
myelin = bwareaopen(darkThings, 1e5);

end


