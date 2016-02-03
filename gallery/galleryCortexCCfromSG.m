function galleryCortexCCfromSG(p, component, outputFile, outputFileSmall)

    % Fix component if necessary
    component(component == 0) = [];
    component = unique(component);

    % Create 3D-spherical structuring element
    rT = 1;
    [x,y,z] = meshgrid(-rT:rT,-rT:rT,-rT:rT);
    se = (x/rT).^2 + (y/rT).^2 + (z/rT).^2 <= 1; 

    % Determine cube (as in p.local) of each ID in component
    load([p.saveFolder 'segToCubeMap.mat']);
    cubeIdx = segToCubeMap(component);
    uniqueCubeIdx = unique(cubeIdx);

    % calculate local isosurfaces in global coordinates 
    for i=1:length(uniqueCubeIdx)
        idx = cubeIdx == uniqueCubeIdx(i);
        theseSegId = component(idx);
        % Display progress
        display(['Processing cube: ' num2str(i) '/' num2str(length(uniqueCubeIdx))]);
        % Load segmentation (global)
        load([p.local(uniqueCubeIdx(i)).saveFolder 'segGlobal.mat']);
        % Restrict to bboxSmall (plus two pixel overlap for continuity of isosurfaces)
        seg = seg(255:end-254,255:end-254,127:end-126);
        % Pad with 2 pixel to include 'end caps'
        seg = padarray(seg, [2 2 2]);
        % Offset for current cube
        zeroOfCube = p.local(uniqueCubeIdx(i)).bboxSmall(:,1)' - [5 5 5];;
        % collect supervoxel
        cube = zeros(size(seg), 'single');
        for k=1:length(theseSegId)
            cube(seg == theseSegId(k)) = 1;
        end
        % Fix border between segments & smooth
        cube = imclose(cube, se);
        cube = smooth3(cube, 'gaussian', 5, 2);
        % Calculate isosurfaces
        issfs{i} = isosurface(cube, .2);
        % Put in nm-space
        issfs{i}.vertices(:,[1 2]) = issfs{i}.vertices(:,[2 1]);
        issfs{i}.vertices = bsxfun(@plus, issfs{i}.vertices, zeroOfCube);
        issfs{i}.vertices = bsxfun(@times, issfs{i}.vertices, [11.24 11.24 28]);
        % Reduce due to memory problems for larger cells
        issfs{i} = reducepatch(issfs{i}, .1);
        issfsSmall{i} = reducepatch(issfs{i}, .1);
    end
    % Make directory if it does not exist
    direc = fileparts(outputFile); 
    if ~exist(direc, 'dir')
        mkdir(direc);
    end
    % Save
    exportSurfaceToAmira(fixIsoStruct(issfs), outputFile);
    exportSurfaceToAmira(fixIsoStruct(issfsSmall), outputFileSmall);
end

function issfsNew = fixIsoStruct(issfs)
    temp = cat(1, issfs{:});
    nrVertices = cellfun(@(x)size(x.vertices,1), issfs);
    vertexOffset = cumsum([0 nrVertices(1:end-1)]);
    for i=1:length(temp)
        temp(i).faces = temp(i).faces + vertexOffset(i);
    end
    issfsNew{1}.vertices = cat(1, temp(:).vertices);
    issfsNew{1}.faces = cat(1, temp(:).faces); 
end

