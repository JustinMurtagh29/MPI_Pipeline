function galleryCortexCCfromSG_new(p, component, outputFile, resolution)

    % Fix component
    component(component == 0) = [];
    component = unique(component);

    % Determine cube (as in p.local) of each ID in component
    load([p.saveFolder 'segToCubeMap.mat']);
    cubeIdx = segToCubeMap(component);
    uniqueCubeIdx = unique(cubeIdx);

    issfs = cell(size(component));
    counter = 1;
    % calculate local isosurfaces in global coordinates 
    for i=1:length(uniqueCubeIdx)
        % Display progress
        display(['Processing cube: ' num2str(i) '/' num2str(length(uniqueCubeIdx))]);
        % Get segId in cubes 
        idx = cubeIdx == uniqueCubeIdx(i);
        theseSegId = unique(component(idx)); 
        % collect isosurfaces of supervoxel
        for k=1:length(theseSegId)
            temp = load([p.local(uniqueCubeIdx(i)).saveFolder '/segments/' num2str(theseSegId(k)) '/' 'isoSurf.mat'], resolution);
            load([p.local(uniqueCubeIdx(i)).saveFolder '/segments/' num2str(theseSegId(k)) '/' 'mask.mat'], 'boxGlobal');
            temp.(resolution).vertices = bsxfun(@plus, temp.(resolution).vertices, double(boxGlobal(:,1)'));
            temp.(resolution).vertices = bsxfun(@times, temp.(resolution).vertices, [11.24 11.24 28]);
            issfs{counter} = temp.(resolution);
            counter = counter + 1;
        end
    end
    % Make directory if it does not exist
    direc = fileparts(outputFile); 
    if ~exist(direc, 'dir')
        mkdir(direc);
    end
    % Save
    exportSurfaceToAmira(issfs', outputFile);
    save(strrep(outputFile, '.issf', '.mat'), 'issfs');
end

