function galleryCortexCCfromSG(p, component, outputFile)

    % Fix component
    component(component == 0) = [];
    component = unique(component);

    % Determine cube (as in p.local) of each ID in component
    load([p.saveFolder 'segToCubeMap.mat']);
    cubeIdx = segToCubeMap(component);
    uniqueCubeIdx = unique(cubeIdx);

    % calculate local isosurfaces in global coordinates 
    for i=1:length(uniqueCubeIdx)
        idx = cubeIdx == uniqueCubeIdx(i);
        theseSegId = unique(component(idx));
        % Display progress
        display(['Processing cube: ' num2str(i) '/' num2str(length(uniqueCubeIdx))]);
        % Load segmentation (global)
        load([p.local(uniqueCubeIdx(i)).saveFolder 'segGlobal.mat']);
        % Initalize variables 
        zeroOfCube = p.local(uniqueCubeIdx(i)).bboxBig';
        % collect supervoxel and calculate isosurfaces
        if ~isempty(theseSegId)
            cube = false(size(seg));
            for k=1:length(theseSegId)
                cube(seg == theseSegId(k)) = true;
            end
            cube = imclose(cube, ones([3,3,3]));
            cube = padarray(cube, [2 2 2]);
            cube = smooth3(cube, 'gaussian', 5, 2);
            issfs{i} = isosurface(cube, .2);
            if ~isempty(issfs{i}.vertices)
                issfs{i} = reducepatch(issfs{i}, .01);
                issfs{i}.vertices(:,[1 2]) = issfs{i}.vertices(:,[2 1]); 			    
                issfs{i}.vertices = bsxfun(@plus, issfs{i}.vertices, zeroOfCube(1,:) - [2 2 2]);
                issfs{i}.vertices = bsxfun(@times, issfs{i}.vertices, [11.24 11.24 28]);
            end
        end
    end
    if exist('issfs', 'var')
        % Remove empty isosurfaces
        idx = false(length(issfs),1);
        for i=1:length(issfs)
            if isempty(issfs{i}) || isempty(issfs{i}.vertices)
                idx(i) = true;
            end
        end
        issfs(idx) = [];
        % Make directory if it does not exist
        direc = fileparts(outputFile); 
        if ~exist(direc, 'dir')
            mkdir(direc);
        end
        % Save
        exportSurfaceToAmira(issfs, outputFile);
        save(strrep(outputFile, '.issf', '.mat'), 'issfs');
        % Reduce patches again
        for i=1:length(issfs)
            issfs{i} = reducepatch(issfs{i}, .01);
        end
        % Remove empty isosurfaces
        idx = false(length(issfs),1);
        for i=1:length(issfs)
            if isempty(issfs{i}) || isempty(issfs{i}.vertices)
                idx(i) = true;
            end
        end
        issfs(idx) = [];
        exportSurfaceToAmira(issfs, strrep(outputFile, 'larger', 'smaller'));
    end
end

