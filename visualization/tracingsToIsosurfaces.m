function tracingsToIsosurfaces(p, source, target, voxelSize, prefix)
    % TRACINGSTOISOSURFACES
    %   a wrapper for Alessandro's isosurface functions. Turn a
    %   folder with nml tracings to a mat file with isosurfaces already corrected for anisotropy.
    %   Intented for use on cluster.
    %  
    %
    % p
    %   Parameters produced by `run configuration.m`
    %
    % prefix
    %   a prefix of the nml files e.g.
    %   're4_RN_CC_S1BFv3__explorational__rkostal__', but a few characters
    %   are enough.
    %   
    %
    % source
    %   path to folder with tracings
    %
    % target
    %   path to output folder
    %
    % voxelSize in nanometers
    %   a row vector e.g. [12 12 30]
    
    % create a job
    isocalc = createJob(CLUSTER_CPU);
    
    % create a task for each nml file
    [~,names] = fileattrib(strcat(source,prefix,'*'));
    
    for i = 1:length(names)
        thisNml = names(1).Name;
        createTask(isocalc, @Visualization.buildIsoSurfaceOfSkel, 1, {p, thisNml});
    end;
    
    % submit
    submit(isocalc);
    wait(isocalc);
   
    
    % correct for anisotropy
    anisotro = fetchOutputs(isocalc);
    
    for i = 1:length(g)
        isotro{i,1}{1,1} = correctAnisotropy(g{i,1}{1,1}, voxelSize);
    end
   
    
    % save in output folder
       
    save(strcat(target,'isosFromTracings.mat'),'isotro')
   


    function outMeshStruct = correctAnisotropy(meshStruct, voxelSize)
        % CORRECTANISOTROPY Scales a mesh according to the voxel size
        % passed as argument. This enables the conversion of pixel
        % numbers to physical units (e.g. nano-metres)
        %
        % Written by
        %   Alessandro Motta <alessandro.motta@brain.mpg.de>
        verts = meshStruct.vertices;
        faces = meshStruct.faces;
    
        % make sure voxel size is a row vector
        voxelSize = reshape(voxelSize, 1, 3);
    
        % multiply indices with direction-dependent
        % measure of nano-metres per pixel
        verts = bsxfun(@times, verts, voxelSize);
    
        % build new mesh
        outMeshStruct = struct;
        outMeshStruct.vertices = verts;
        outMeshStruct.faces = faces;
        end
    
        
end