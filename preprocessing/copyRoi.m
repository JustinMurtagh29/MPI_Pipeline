function raw = copyRoi()
    % Copy part of 2012_09_28_ex145_07x2_new2 dataset availible on webKnossos to a new dataset to have the ROI for segmentation aligned to KNOSSOS cubes 

    % Define source of copy process
    chosenOne.root = '/gaba/u/kboerg/st07x2_new2/color/1/';
    chosenOne.prefix = '2012-09-28_ex145_07x2_new2_mag1'; 
    chosenOne.bbox = [1195, 1515, 115, 6690, 9945, 3420];
    % Convert to Matlab format (and add +1 offset between coordinate systems)
    chosenOne.bbox = Util.convertWebknossosToMatlabBbox(chosenOne.bbox);
    
    % Define destination of copy process
    destination.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016/color/1/';
    destination.prefix = '2012-09-28_ex145_07x2_ROI2016_mag1';
    destination.firstCoord = 129 - [25 25 10];

    % Read from one, write to wKcubes
    raw = readKnossosRoi(chosenOne.root, chosenOne.prefix, chosenOne.bbox);
    writeKnossosRoi(destination.root, destination.prefix, destination.firstCoord, raw);

    % Saving raw data to HDF 5 file for faster loading in case we need multiple tries for contrast normalization
    % Saving duplicates in memory not feasible
    h5create('/gaba/scratch/mberning/tempData07x2.h5', '/raw', size(raw), 'Datatype', 'uint8');
    h5write('/gaba/scratch/mberning/tempData07x2.h5', '/raw', raw);

end

