function bloodVesselVisualization(bloodVessels, bbox)
    % Visualize blood vessels for 07x2 from dataset preprocessing

    z = bbox(3,1):50:bbox(3,2);
    for i=1:length(z)-1
        display(['-----' num2str(i) '/' num2str(length(z)-1) '-----']);
        display('Smoothing blood vessels');
        tic;
        smoothedPart = smooth3(bloodVessels(:,:,z(i):z(i+1)+10), 'gaussian', 5, 2);
        smoothedPart = padarray(smoothedPart, [2 2 2], 0);
        toc;
        display('Isosurface calculation');
        tic;
        issfs{i} = isosurface(smoothedPart, .1);
        toc;
        display('Simplyfing isosurfaces');
        tic;
        issfs{i} = reducepatch(issfs{i}, .05);
        toc;
        display('Scaling isosurfaces');
        tic;
        % Position isosurfaces in nm at right position to be consitent with normal dataset coordinates 
        issfs{i}.vertices = bsxfun(@plus, issfs{i}.vertices, [bbox(1:2,1); z(i)]' - [2 2 2]);
        issfs{i}.vertices = bsxfun(@times, issfs{i}.vertices, [11.24 11.24 28]);
        toc;
    end

    save('/zdata/manuel/results/bloodVesselIsosurface.mat', 'issfs');
    exportSurfaceToAmira(issfs, '/zdata/manuel/results/bloodVesselsIsosurface.issf');

end

