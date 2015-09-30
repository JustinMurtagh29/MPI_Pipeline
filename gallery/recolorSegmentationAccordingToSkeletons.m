function recolorSegmentationAccordingToSkeletons()
    
    cubes = [8 4 2; 29 38 41];
    skelPath = '/u/mberning/data/retina/retinaN2skeletons/allNice/';
    segParam.root = '/gaba/u/mberning/backup/nfs_mberning/20140310backup/mag1/'; 
    segParam.prefix = '100527_k0563_seg'; 
    segOutParam.root = '/gaba/u/mberning/results/retinaRecoloring/';
    segOutParam.prefix = '100527_k0563_seg';

    % Load skeletons, extract nodes
    files = dir([skelPath '*.nml']);
    for i=1:length(files)
        % Use in this way to supress annoying output
        [~, skel_data(i)] = evalc('parseNml([skelPath files(i).name])');
        nodes{i} = skel_data{i}.nodes(:,1:3);
    end

    % Start one process for each cube (on cluster)
    cubesProcessed = 0;
    for x=cubes(1,1):cubes(2,1)
        for y=cubes(1,2):cubes(2,2)
            for z=cubes(1,3):cubes(2,3)
                cubesProcessed = cubesProcessed + 1;
                inputCell{cubesProcessed} = {segParam, segOutParam, nodes, [x y z]};
            end
        end
    end
    functionH = @recolorCubeAccordingToNodes;
    startCPU(functionH, inputCell, 'retina recoloring');
    save([segOutParam.root 'skeletonToColor.mat'], 'files', 'skel_data');

end

