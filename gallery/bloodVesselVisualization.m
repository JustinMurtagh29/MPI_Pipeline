function bloodVesselVisualization(bloodVessels, bbox)
% Visualize blood vessels for 07x2 from dataset preprocessing

display('Smoothing blood vessels');
tic;
bloodVessels = smooth3(bloodVessels, 'gaussian', 6, 3);
toc;
display('Isosurface calculation');
tic;
issfs = isosurface(bloodVessels, .1);
toc;
display('Simplyfing isosurfaces');
tic;
issfs = reducepatch(issfs, .01);
toc;
display('Scaling and saving isosurfaces');
tic;
issfs.vertices = bsxfun(@plus, issfs.vertices, bbox(:,1)' + [1 1 1]);
issfs.vertices = bsxfun(@times, issfs.vertices, [11.24 11.24 28]);
save('/zdata/manuel/results/bloodVesselIsosurface.mat', 'issfs');
exportSurfaceToAmira(issfs, '/zdata/manuel/results/bloodVesselsIsosurface.issf');
toc;

end
