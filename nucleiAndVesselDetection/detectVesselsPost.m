function vesselsPost = detectVesselsPost(vessels)
tic;
display('Postprocessing');
% Smooth and connect surface (by closing with spherical (in nm) strucuturing element)
[x,y,z] = meshgrid(-12:12,-12:12,-5:5);
se = (x/12).^2 + (y/12).^2 + (z/5).^2 <= 1;
vesselsPost = imclose(vessels, se);
% Remove very small objects (smaller than 100*10*10 voxel
vesselsPost = bwareaopen(vesselsPost, 10000);
toc;
